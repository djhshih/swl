import os

import swl.ir.node as ir
from swl.semantic.task.type import Param, TaskSignature
from swl.semantic.wf.check import Checker
from swl.syntax.wf import builtins, node as wf_node


def _builtin_map_signature():
    return TaskSignature({'f': Param('f', None), 'xs': Param('xs', None)}, {}, {})


def _builtin_map_by_signature():
    return TaskSignature({'f': Param('f', None), 'key': Param('key', None), 'xs': Param('xs', None)}, {}, {})


def _normalize_twice(lowerer, node):
    return lowerer.normalize(lowerer.normalize(node))


def _literal_or_none(expr):
    return expr.value if expr.type == wf_node.NodeType.str else None


def _extend_signature_inputs(inputs, available, sig):
    for name, param in sig.inputs.items():
        if name not in available:
            inputs[name] = param
    available.update(sig.inputs.keys())
    available.update(sig.outputs.keys())


def _lower_block_items(lowerer, items, env, imports):
    bindings = []
    local_env = dict(env)
    for item in items[:-1]:
        if item.type != wf_node.NodeType.bind:
            continue
        value = lowerer.lower_binding(item, local_env, imports)
        local_env = dict(local_env)
        if item.id.name in imports:
            local_env[item.id.name] = value
            continue
        var = ir.Variable(lowerer._alloc_var_id(), item.id.name, value)
        local_env[item.id.name] = ir.Ref(var.id, var.name)
        bindings.append(var)
    return bindings, local_env, lowerer.lower_expr(items[-1], local_env, imports)


class Lowerer:
    def __init__(self, files=None):
        self.checker = Checker(files=files)
        self.workflow_cache = {}
        self.function_cache = {}
        self.next_var_id = 1
        self.next_generated_id = 1
        self.signature = None

    def lower_file(self, path: str):
        result = self.checker.load(path)
        if result.errors:
            raise ValueError('\n'.join(result.errors))
        self.signature = result.signature
        return self.lower_tree(result.tree, result.imports, result.signature)

    def lower_tree(self, tree, imports, signature=None):
        old_sig = self.signature
        self.signature = signature
        try:
            return self._lower_tree_impl(tree, imports, signature)
        finally:
            self.signature = old_sig

    def _lower_tree_impl(self, tree, imports, signature=None):
        env = {}
        if tree.type != wf_node.NodeType.block:
            return self.normalize(self._ensure_block_root(self.lower_expr(tree, env, imports), signature))

        bindings, env, result = _lower_block_items(self, tree.body, env, imports)
        if isinstance(result, ir.Lambda) and signature is not None:
            result = ir.Lambda(result.param, result.body, signature, getattr(result, 'is_batch', False))
        if not bindings:
            return _normalize_twice(self, result)
        return _normalize_twice(self, ir.Block(bindings, result))

    def lower_binding(self, expr, env, imports):
        if expr.type == wf_node.NodeType.bind and expr.id.name in imports:
            imported = imports[expr.id.name]
            return self._function_from_import(expr.id.name, imported)
        return self.lower_expr(expr.value, env, imports)

    def _lower_name(self, name, env, imports):
        if name in env:
            return env[name]
        if name in imports:
            return self._function_from_import(name, imports[name])
        if self.signature is not None and name in self.signature.inputs:
            param = self.signature.inputs[name]
            typ = param.type.value if param.type is not None else None
            return ir.Input(name, typ, param.desc)
        raise ValueError(f'Undefined variable: {name}')

    def lower_expr(self, expr, env, imports):
        if expr.type == wf_node.NodeType.id:
            return self._lower_name(expr.name, env, imports)

        if expr.type in (wf_node.NodeType.num, wf_node.NodeType.str):
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
            import_path = builtins.match_import(expr)
            if import_path is not None:
                return self._lower_inline_import(import_path)
            builtin = self._match_builtin(expr)
            if builtin == 'map_apply':
                function_expr, arg_expr = builtins.match_map(expr)
                return ir.Map(
                    self.lower_expr(function_expr, env, imports),
                    self.lower_expr(arg_expr, env, imports),
                )
            if builtin == 'map_by_apply':
                key_expr, function_expr, arg_expr = builtins.match_map_by(expr)
                return ir.Map(
                    self.lower_expr(function_expr, env, imports),
                    self.lower_expr(arg_expr, env, imports),
                    key=_literal_or_none(key_expr),
                )
            if builtin == 'map_partial':
                return ir.Apply(
                    ir.Function('map', 'builtin', _builtin_map_signature()),
                    self.lower_expr(expr.arg, env, imports),
                )
            if builtin == 'map_by_partial':
                return ir.Apply(
                    ir.Function('map_by', 'builtin', _builtin_map_by_signature()),
                    self.lower_expr(expr.arg, env, imports),
                )
            satisfied = self.checker._apply_satisfied.get(id(expr), set())
            return ir.Apply(
                self.lower_expr(expr.fun, env, imports),
                self.lower_expr(expr.arg, env, imports),
                satisfied=satisfied,
            )

        if expr.type == wf_node.NodeType.fun:
            local_env = dict(env)
            local_env[expr.param.name] = ir.Name(expr.param.name)
            body = self.lower_expr(expr.body, local_env, imports)
            return ir.Lambda(expr.param.name, self._ensure_block(body))

        if expr.type == wf_node.NodeType.block:
            bindings, _, result = _lower_block_items(self, expr.body, env, imports)
            if not bindings:
                return result
            return ir.Block(bindings, result)

        if expr.type == wf_node.NodeType.chain:
            items = self._lower_chain_items(expr, env, imports)
            return self._chain_to_lambda_block(items, None)

        if expr.type == wf_node.NodeType.bind:
            value = self.lower_binding(expr, env, imports)
            return ir.Variable(self._alloc_var_id(), expr.id.name, value)

        raise ValueError(f'Unhandled AST node type during lowering: {expr.type}, expression: {expr}')


    def _lower_inline_import(self, path):
        from swl.semantic.wf.imports import load_import
        base_dir = os.path.dirname(self.checker._loading[-1]) if self.checker._loading else '.'
        full_path = os.path.abspath(os.path.join(base_dir, path))
        stem = os.path.splitext(os.path.basename(full_path))[0]
        imported = load_import(self.checker, stem, full_path)
        return self._function_from_import(stem, imported)

    def _function_from_import(self, name, imported):
        if name in self.function_cache:
            return self.function_cache[name]
        if imported.kind == 'workflow':
            body = self._cached_workflow_body(imported.path)
            function = ir.Function(name, imported.kind, imported.signature, imported.path, body, getattr(imported.check, 'is_batch', False))
        else:
            function = ir.Function(name, imported.kind, imported.signature, imported.path, None, False)
        self.function_cache[name] = function
        return function

    def _cached_workflow_body(self, path):
        if path in self.workflow_cache:
            return self.workflow_cache[path]
        result = self.checker.load(path)
        body = self.lower_tree(result.tree, result.imports, result.signature)
        self.workflow_cache[path] = body
        return body

    def _match_builtin(self, expr):
        if expr.type != wf_node.NodeType.apply:
            return None
        if builtins.match_map(expr) is not None:
            return 'map_apply'
        if builtins.match_map_by(expr) is not None:
            return 'map_by_apply'
        if expr.fun.type == wf_node.NodeType.id and expr.fun.name == 'map':
            return 'map_partial'
        if expr.fun.type == wf_node.NodeType.id and expr.fun.name == 'map_by':
            return 'map_by_partial'
        return None

    def _lower_chain_items(self, expr, env, imports):
        if expr.type != wf_node.NodeType.chain:
            return [self.lower_expr(expr, env, imports)]
        return self._lower_chain_items(expr.left, env, imports) + self._lower_chain_items(expr.right, env, imports)

    def normalize(self, node):
        if isinstance(node, ir.Lambda):
            body = self.normalize(node.body)
            staged = self._normalize_block_lambda(node.param, body, node.signature)
            if staged is not None:
                return staged
            return ir.Lambda(node.param, body, node.signature)

        if isinstance(node, ir.Block):
            raw_bindings = [ir.Variable(bind.id, bind.name, self.normalize(bind.value)) for bind in node.bindings]
            binding_map = {bind.id: bind.value for bind in raw_bindings}
            bindings = [ir.Variable(bind.id, bind.name, self._normalize_with_bindings(bind.value, binding_map)) for bind in raw_bindings]
            binding_map = {bind.id: bind.value for bind in bindings}
            result = self._normalize_with_bindings(node.result, binding_map)
            return ir.Block(bindings, result)

        if isinstance(node, ir.Variable):
            return ir.Variable(node.id, node.name, self.normalize(node.value))

        if isinstance(node, ir.Apply):
            return ir.Apply(self.normalize(node.function), self.normalize(node.arg), node.signature)

        if isinstance(node, ir.Map):
            fn = self.normalize(node.function)
            arg = self.normalize(node.arg)
            fn = self._materialize_mappable_callable(fn)
            return ir.Map(fn, arg, node.key)

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
        if not isinstance(node, ir.Update):
            return [node]
        return self._flatten_update_terms(node.left) + self._flatten_update_terms(node.right)

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
            _extend_signature_inputs(inputs, available, sig)
        outputs = {}
        run = {}
        for sig in signatures:
            outputs.update(sig.outputs)
            run.update(sig.run)
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

    def _ensure_block(self, node):
        return node if isinstance(node, ir.Block) else ir.Block([], node)

    def _ensure_block_root(self, node, signature):
        if isinstance(node, ir.Lambda):
            sig = signature if signature is not None else node.signature
            return ir.Lambda(node.param, self._ensure_block(node.body), sig)
        return node

    def _normalize_with_bindings(self, node, binding_map):
        if isinstance(node, ir.Map):
            fn = self._normalize_with_bindings(node.function, binding_map)
            if isinstance(fn, ir.Ref) and fn.id in binding_map:
                fn = binding_map[fn.id]
            fn = self._materialize_mappable_callable(fn)
            arg = self._normalize_with_bindings(node.arg, binding_map)
            return ir.Map(fn, arg, node.key)
        if isinstance(node, ir.Apply):
            return ir.Apply(self._normalize_with_bindings(node.function, binding_map), self._normalize_with_bindings(node.arg, binding_map), node.signature)
        if isinstance(node, ir.Record):
            return ir.Record({name: self._normalize_with_bindings(value, binding_map) for name, value in node.fields.items()})
        if isinstance(node, ir.Field):
            return ir.Field(self._normalize_with_bindings(node.record, binding_map), node.name)
        if isinstance(node, ir.Update):
            return ir.Update(self._normalize_with_bindings(node.left, binding_map), self._normalize_with_bindings(node.right, binding_map))
        if isinstance(node, ir.Block):
            inner = {bind.id: bind.value for bind in node.bindings}
            merged = dict(binding_map)
            merged.update(inner)
            return ir.Block(node.bindings, self._normalize_with_bindings(node.result, merged))
        return node

    def _materialize_mappable_callable(self, node):
        if isinstance(node, ir.Function):
            return node
        if isinstance(node, ir.Apply):
            function = self._materialize_mappable_callable(node.function)
            if isinstance(function, ir.Function):
                helper = self._generated_callable_from_apply(function, node.arg)
                if helper is not None:
                    return helper
            return ir.Apply(function, node.arg, node.signature)
        if isinstance(node, ir.Lambda):
            helper = self._generated_callable_from_lambda(node)
            if helper is not None:
                return helper
        return node

    def _generated_lambda_signature(self, param, outputs):
        return TaskSignature({param: Param(param, None)}, {out: Param(out, None) for out in outputs}, {})

    def _generated_callable_from_lambda(self, node):
        if not isinstance(node, ir.Lambda):
            return None
        helper = ir.Lambda(node.param, self.normalize(node.body), node.signature)
        outputs = self._infer_lambda_output_names(helper)
        if outputs is None:
            return None
        name = f'map_lambda_{self.next_generated_id}'
        self.next_generated_id += 1
        signature = self._generated_lambda_signature(node.param, outputs)
        return ir.Function(
            name,
            'workflow',
            signature,
            f'<generated:{name}>',
            body=ir.Lambda(node.param, helper.body, signature),
        )

    def _generated_callable_from_apply(self, function, arg):
        if not isinstance(function, ir.Function):
            return None
        if not isinstance(arg, ir.Record):
            return None
        bound_names = set(arg.fields.keys())
        remaining = [name for name in function.signature.inputs.keys() if name not in bound_names]
        if not remaining:
            return None
        from swl.semantic.task.type import Param, TaskSignature
        name = f'map_partial_{self.next_generated_id}'
        self.next_generated_id += 1
        signature = TaskSignature({'x': Param('x', None)}, dict(function.signature.outputs), {})
        return ir.Function(
            name,
            'workflow',
            signature,
            f'<generated:{name}>',
            body=ir.Lambda('x', ir.Block([], ir.Apply(function, ir.Update(ir.Record({name: ir.Field(ir.Name('x'), name) for name in remaining}), arg))), signature),
        )

    def _infer_lambda_output_names(self, node):
        body = node.body.result if isinstance(node.body, ir.Block) else None
        if isinstance(body, ir.Record):
            return list(body.fields.keys())
        if isinstance(body, ir.Apply) and isinstance(body.function, ir.Function):
            return list(body.function.signature.outputs.keys())
        return None

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
