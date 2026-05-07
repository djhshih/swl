import os

from swl.ir import node as ir
from swl.semantic.wf.check import Checker
from swl.syntax.wf import node as wf_node
from swl.syntax.wf.parser import Parser as WfParser


class Lowerer:
    def __init__(self):
        self.checker = Checker()

    def lower_file(self, path: str):
        result = self.checker.load(path)
        return self.lower_tree(result.tree, result.imports, result.signature)

    def lower_tree(self, tree, imports, signature=None):
        env = {}
        if tree.type != wf_node.NodeType.block:
            return self.lower_expr(tree, env, imports)

        bindings = []
        exprs = tree.body
        for expr in exprs[:-1]:
            if expr.type == wf_node.NodeType.bind:
                value = self.lower_binding(expr, env, imports)
                bindings.append(ir.Bind(expr.id.name, value))
                env = dict(env)
                env[expr.id.name] = value
        result = self.lower_expr(exprs[-1], env, imports)
        return ir.Block(bindings, result)

    def lower_expr(self, expr, env, imports):
        if expr.type == wf_node.NodeType.id:
            if expr.name in env:
                return env[expr.name]
            if expr.name in imports:
                imported = imports[expr.name]
                return ir.Import(expr.name, imported.kind, imported.path, imported.signature)
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
            return ir.Lambda(
                expr.param.name,
                self.lower_expr(expr.body, local_env, imports),
            )

        if expr.type == wf_node.NodeType.block:
            bindings = []
            local_env = dict(env)
            for item in expr.body[:-1]:
                if item.type == wf_node.NodeType.bind:
                    value = self.lower_binding(item, local_env, imports)
                    bindings.append(ir.Bind(item.id.name, value))
                    local_env = dict(local_env)
                    local_env[item.id.name] = value
            result = self.lower_expr(expr.body[-1], local_env, imports)
            return ir.Block(bindings, result)

        if expr.type == wf_node.NodeType.chain:
            return ir.Chain(self._lower_chain_items(expr, env, imports))

        if expr.type == wf_node.NodeType.bind:
            return ir.Bind(expr.id.name, self.lower_expr(expr.value, env, imports))

        return ir.Unknown()

    def lower_binding(self, expr, env, imports):
        if expr.type == wf_node.NodeType.bind and expr.id.name in imports:
            imported = imports[expr.id.name]
            return ir.Import(expr.id.name, imported.kind, imported.path, imported.signature)
        return self.lower_expr(expr.value, env, imports)

    def _lower_chain_items(self, expr, env, imports):
        if expr.type == wf_node.NodeType.chain:
            return self._lower_chain_items(expr.left, env, imports) + self._lower_chain_items(expr.right, env, imports)
        return [self.lower_expr(expr, env, imports)]


def lower_file(path: str):
    return Lowerer().lower_file(path)


def parse_and_lower(src: str, base_dir: str = '.'):
    tree = WfParser().parse(src)
    checker = Checker()
    fake_path = os.path.join(base_dir, '__input__.swl')
    imports = checker._load_imports(tree, os.path.dirname(fake_path))
    return Lowerer().lower_tree(tree, imports)
