import os

from swl.ir import node as ir
from swl.semantic.wf.check import Checker
from swl.syntax.wf import node as wf_node
from swl.syntax.wf.parser import Parser as WfParser


class Lowerer:
    def __init__(self, files=None):
        self.checker = Checker(files=files)
        self.workflow_cache = {}
        self.function_cache = {}
        self.next_var_id = 1

    def lower_file(self, path: str):
        result = self.checker.load(path)
        if result.errors:
            raise ValueError('\n'.join(result.errors))
        return self.lower_tree(result.tree, result.imports, result.signature)

    def lower_tree(self, tree, imports, signature=None):
        env = {}
        if tree.type != wf_node.NodeType.block:
            return self.normalize(self._ensure_block_root(self.lower_expr(tree, env, imports), signature))

        bindings = []
        exprs = tree.body
        for expr in exprs[:-1]:
            if expr.type == wf_node.NodeType.bind:
                value = self.lower_binding(expr, env, imports)
                env = dict(env)
                if expr.id.name in imports:
                    env[expr.id.name] = value
                    continue
                var = ir.Variable(self._alloc_var_id(), expr.id.name, value)
                env[expr.id.name] = ir.Ref(var.id, var.name)
                bindings.append(var)
        result = self.lower_expr(exprs[-1], env, imports)
        if isinstance(result, ir.Lambda):
            if signature is not None:
                result = ir.Lambda(result.param, result.body, signature)
            if not bindings:
                return self.normalize(result)
        if not bindings:
            return self.normalize(result)
        return self.normalize(ir.Block(bindings, result))

    def lower_binding(self, expr, env, imports):
        if expr.type == wf_node.NodeType.bind and expr.id.name in imports:
            imported = imports[expr.id.name]
            return self._function_from_import(expr.id.name, imported)
        return self.lower_expr(expr.value, env, imports)

    def lower_expr(self, expr, env, imports):
        if expr.type == wf_node.NodeType.id:
            if expr.name in env:
                return env[expr.name]
            if expr.name in imports:
                return self._function_from_import(expr.name, imports[expr.name])
            return ir.Name(expr.name)

        if expr.type == wf_node.NodeType.num:
            return ir.Literal(expr.value)

        if expr.type == wf_node.NodeType.str:
            return ir.Literal(expr.value)

        if expr.type == wf_node.NodeType.rec:
            fields = {}
            for name, value in expr.value.items():
                fields[name] = self.lower_expr(value, env, imports)
            return ir.Record(fields)

        if expr.type == wf_node.NodeType.get:
            return ir.Field(self.lower_expr(expr.rec, env, imports), expr.member.name)

        if expr.type == wf_node.NodeType.update:
            return ir.Update(
                self.lower_expr(expr.left, env, imports),
                self.lower_expr(expr.right, env, imports),
            )

        if expr.type == wf_node.NodeType.apply:
            return ir.Apply(
                self.lower_expr(expr.fun, env, imports),
                self.lower_expr(expr.arg, env, imports),
            )

        if expr.type == wf_node.NodeType.fun:
            local_env = dict(env)
            local_env[expr.param.name] = ir.Name(expr.param.name)
            body = self.lower_expr(expr.body, local_env, imports)
            if isinstance(body, ir.Block):
                return ir.Lambda(expr.param.name, body)
            return ir.Lambda(expr.param.name, ir.Block([], body))

        if expr.type == wf_node.NodeType.block:
            bindings = []
            local_env = dict(env)
            for item in expr.body[:-1]:
                if item.type == wf_node.NodeType.bind:
                    value = self.lower_binding(item, local_env, imports)
                    local_env = dict(local_env)
                    if item.id.name in imports:
                        local_env[item.id.name] = value
                        continue
                    var = ir.Variable(self._alloc_var_id(), item.id.name, value)
                    local_env[item.id.name] = ir.Ref(var.id, var.name)
                    bindings.append(var)
            result = self.lower_expr(expr.body[-1], local_env, imports)
            if not bindings:
                return result
            return ir.Block(bindings, result)

        if expr.type == wf_node.NodeType.chain:
            items = self._lower_chain_items(expr, env, imports)
            return self._chain_to_lambda_block(items, None)

        if expr.type == wf_node.NodeType.bind:
            value = self.lower_binding(expr, env, imports)
            return ir.Variable(self._alloc_var_id(), expr.id.name, value)

        return ir.Unknown()

    def _function_from_import(self, name, imported):
        if name in self.function_cache:
            return self.function_cache[name]
        if imported.kind == 'workflow':
            body = self._cached_workflow_body(imported.path)
            function = ir.Function(name, imported.kind, imported.signature, imported.path, body)
        else:
            function = ir.Function(name, imported.kind, imported.signature, imported.path, None)
        self.function_cache[name] = function
        return function

    def _cached_workflow_body(self, path):
        if path in self.workflow_cache:
            return self.workflow_cache[path]
        result = self.checker.load(path)
        body = self.lower_tree(result.tree, result.imports, result.signature)
        self.workflow_cache[path] = body
        return body

    def _lower_chain_items(self, expr, env, imports):
        if expr.type == wf_node.NodeType.chain:
            return self._lower_chain_items(expr.left, env, imports) + self._lower_chain_items(expr.right, env, imports)
        return [self.lower_expr(expr, env, imports)]

    def normalize(self, node):
        if isinstance(node, ir.Lambda):
            body = self.normalize(node.body)
            staged = self._normalize_block_lambda(node.param, body, node.signature)
            if staged is not None:
                return staged
            return ir.Lambda(node.param, body, node.signature)

        if isinstance(node, ir.Block):
            bindings = [ir.Variable(bind.id, bind.name, self.normalize(bind.value)) for bind in node.bindings]
            result = self.normalize(node.result)
            return ir.Block(bindings, result)

        if isinstance(node, ir.Variable):
            return ir.Variable(node.id, node.name, self.normalize(node.value))

        if isinstance(node, ir.Apply):
            return ir.Apply(self.normalize(node.function), self.normalize(node.arg), node.signature)

        if isinstance(node, ir.Ref):
            return node

        if isinstance(node, ir.Record):
            return ir.Record({name: self.normalize(value) for name, value in node.fields.items()})

        if isinstance(node, ir.Field):
            return ir.Field(self.normalize(node.record), node.name)

        if isinstance(node, ir.Update):
            return self._normalize_update(self.normalize(node.left), self.normalize(node.right))

        return node

    def _normalize_update(self, left, right):
        terms = self._flatten_update_terms(left) + self._flatten_update_terms(right)
        return self._merge_terms(terms)

    def _flatten_update_terms(self, node):
        if isinstance(node, ir.Update):
            return self._flatten_update_terms(node.left) + self._flatten_update_terms(node.right)
        return [node]

    def _merge_terms(self, terms):
        if not terms:
            return ir.Record({})
        result = terms[0]
        for term in terms[1:]:
            result = ir.Update(result, term)
        return result

    def _normalize_block_lambda(self, param, body, signature):
        if not isinstance(body, ir.Block):
            return None
        return ir.Lambda(param, body, signature)

    def _compose_signature(self, items, fallback):
        if fallback is not None:
            return fallback
        signatures = []
        for item in items:
            sig = getattr(item, 'signature', None)
            if sig is None:
                return None
            signatures.append(sig)
        if not signatures:
            return None
        inputs = dict(signatures[0].inputs)
        available = set(signatures[0].inputs.keys()).union(signatures[0].outputs.keys())
        for sig in signatures[1:]:
            for name, param in sig.inputs.items():
                if name not in available:
                    inputs[name] = param
            available.update(sig.inputs.keys())
            available.update(sig.outputs.keys())
        outputs = {}
        run = {}
        for sig in signatures:
            outputs.update(sig.outputs)
            run.update(sig.run)
        from swl.semantic.task.type import TaskSignature
        return TaskSignature(inputs, outputs, run)

    def _chain_to_lambda_block(self, items, signature):
        param = '_input'
        root = ir.Name(param)
        bindings = []
        available = [root]
        for i, item in enumerate(items):
            name = f'_s{i + 1}'
            arg = self._merge_terms(list(available))
            var = ir.Variable(self._alloc_var_id(), name, ir.Apply(item, arg))
            bindings.append(var)
            available.append(ir.Ref(var.id, var.name))
        result = self._merge_terms([ir.Ref(bind.id, bind.name) for bind in bindings]) if bindings else ir.Record({})
        return ir.Lambda(param, ir.Block(bindings, result), self._compose_signature(items, signature))

    def _ensure_block_root(self, node, signature):
        if isinstance(node, ir.Lambda):
            body = node.body if isinstance(node.body, ir.Block) else ir.Block([], node.body)
            sig = signature if signature is not None else node.signature
            return ir.Lambda(node.param, body, sig)
        return node

    def _alloc_var_id(self):
        current = self.next_var_id
        self.next_var_id += 1
        return current


def lower_file(path: str, files=None):
    return Lowerer(files=files).lower_file(path)


def parse_and_lower(src: str, base_dir: str = '.', files=None):
    fake_path = os.path.abspath(os.path.join(base_dir, '__input__.swl'))
    checker = Checker(files=files)
    result = checker.load_content(src, fake_path)
    if result.errors:
        raise ValueError('\n'.join(result.errors))
    return Lowerer(files=files).lower_tree(result.tree, result.imports, result.signature)
