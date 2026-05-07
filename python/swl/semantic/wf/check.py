import os

from swl.semantic.task.type import Param, TaskSignature, TypeChecker, signature_from_task
from swl.syntax.task.parser import Parser as TaskParser
from swl.syntax.wf import node as wf_node
from swl.syntax.wf.parser import Parser as WfParser


class Import:
    def __init__(self, name: str, path: str, signature: TaskSignature, kind: str):
        self.name = name
        self.path = path
        self.signature = signature
        self.kind = kind


class OpenRecord:
    def __init__(self, fields=None):
        self.fields = set(fields or [])


class ClosedRecord:
    def __init__(self, fields=None):
        self.fields = set(fields or [])


class TaskResult:
    def __init__(self, task_name: str, signature: TaskSignature):
        self.task_name = task_name
        self.signature = signature


class UnknownValue:
    pass


class WorkflowCheck:
    def __init__(self, tree, imports, chain_errors, inferred_inputs, issues, signature=None):
        self.tree = tree
        self.imports = imports
        self.chain_errors = chain_errors
        self.inferred_inputs = inferred_inputs
        self.issues = issues
        self.signature = signature


class Checker:
    def __init__(self):
        self._loading = []

    def load(self, path: str) -> WorkflowCheck:
        full_path = os.path.abspath(path)
        if full_path in self._loading:
            raise ValueError(f'Circular workflow import: {full_path}')
        self._loading.append(full_path)
        try:
            with open(full_path, 'r') as f:
                src = f.read()
            tree = WfParser().parse(src)
            imports = self._load_imports(tree, os.path.dirname(full_path))
            checker = TypeChecker()
            for imported in imports.values():
                checker.add_task(imported.name, imported.signature)
            chain_errors = self._check_chains(tree, checker)
            inferred_inputs, issues = self._infer_inputs(tree, imports)
            signature = self._build_workflow_signature(tree, imports, inferred_inputs, issues)
            return WorkflowCheck(tree, imports, chain_errors, inferred_inputs, issues, signature)
        finally:
            self._loading.pop()

    def _load_imports(self, tree, base_dir: str):
        imports = {}
        if tree.type != wf_node.NodeType.block:
            return imports
        for expr in tree.body:
            if expr.type != wf_node.NodeType.bind:
                continue
            path = self._match_import(expr.value)
            if path is None:
                continue
            full_path = os.path.join(base_dir, path)
            imports[expr.id.name] = self._load_import(expr.id.name, full_path)
        return imports

    def _load_import(self, name: str, path: str) -> Import:
        if path.endswith('.sh'):
            with open(path) as f:
                task = TaskParser().parse(f.read())
            return Import(name, path, signature_from_task(task), 'task')
        if path.endswith('.swl'):
            check = self.load(path)
            if check.signature is None:
                raise ValueError(f'Imported workflow does not produce a signature: {path}')
            return Import(name, path, check.signature, 'workflow')
        raise ValueError(f'Unrecognized import path: {path}')

    def _match_import(self, expr):
        if expr.type != wf_node.NodeType.apply:
            return None
        if expr.fun.type != wf_node.NodeType.id or expr.fun.name != 'import':
            return None
        if expr.arg.type != wf_node.NodeType.str:
            return None
        return expr.arg.value

    def _check_chains(self, tree, checker: TypeChecker):
        errors = []
        self._walk_chains(tree, checker, errors)
        return errors

    def _walk_chains(self, expr, checker: TypeChecker, errors):
        if expr.type == wf_node.NodeType.chain:
            left_name = self._task_name(expr.left)
            right_name = self._task_name(expr.right)
            if left_name and right_name:
                errors.extend(checker.check_chain(left_name, right_name))
        for child in self._children(expr):
            self._walk_chains(child, checker, errors)

    def _infer_inputs(self, tree, imports):
        if tree.type != wf_node.NodeType.block or not tree.body:
            return set(), []

        final = tree.body[-1]
        issues = []
        demanded = set()

        if final.type == wf_node.NodeType.fun:
            env = {final.param.name: OpenRecord()}
            self._eval_function_body(final, imports, env, demanded, issues)
            return demanded, issues

        env = {}
        for expr in tree.body[:-1]:
            if expr.type == wf_node.NodeType.bind:
                env[expr.id.name] = self._eval_expr(expr.value, imports, env, demanded, issues)
        self._eval_expr(final, imports, env, demanded, issues)
        return demanded, issues

    def _build_workflow_signature(self, tree, imports, inferred_inputs, issues):
        if tree.type != wf_node.NodeType.block or not tree.body:
            return None

        final = tree.body[-1]
        if final.type == wf_node.NodeType.fun:
            env = {final.param.name: OpenRecord(inferred_inputs)}
            result = self._eval_function_body(final, imports, env, set(), issues)
        else:
            env = {}
            for expr in tree.body[:-1]:
                if expr.type == wf_node.NodeType.bind:
                    env[expr.id.name] = self._eval_expr(expr.value, imports, env, set(), issues)
            result = self._eval_expr(final, imports, env, set(), issues)

        inputs = {
            name: Param(name, None)
            for name in sorted(inferred_inputs)
        }
        outputs = self._signature_outputs(result)
        return TaskSignature(inputs, outputs, {})

    def _signature_outputs(self, value):
        if isinstance(value, TaskResult):
            return dict(value.signature.outputs)
        if isinstance(value, ClosedRecord):
            return {name: Param(name, None) for name in sorted(value.fields)}
        if isinstance(value, OpenRecord):
            return {name: Param(name, None) for name in sorted(value.fields)}
        return {}

    def _eval_function_body(self, fn, imports, env, demanded, issues):
        local_env = dict(env)
        if fn.body.type != wf_node.NodeType.block:
            return self._eval_expr(fn.body, imports, local_env, demanded, issues)

        result = UnknownValue()
        for expr in fn.body.body:
            if expr.type == wf_node.NodeType.bind:
                value = self._eval_expr(expr.value, imports, local_env, demanded, issues)
                local_env[expr.id.name] = value
                result = value
            else:
                result = self._eval_expr(expr, imports, local_env, demanded, issues)
        return result

    def _eval_expr(self, expr, imports, env, demanded, issues):
        if expr.type == wf_node.NodeType.id:
            return env.get(expr.name, UnknownValue())

        if expr.type == wf_node.NodeType.rec:
            fields = set(expr.value.keys())
            for value_expr in expr.value.values():
                self._eval_expr(value_expr, imports, env, demanded, issues)
            return ClosedRecord(fields)

        if expr.type == wf_node.NodeType.get:
            rec = self._eval_expr(expr.rec, imports, env, demanded, issues)
            field = expr.member.name
            if isinstance(rec, OpenRecord):
                demanded.add(field)
                rec.fields.add(field)
                return UnknownValue()
            if isinstance(rec, ClosedRecord):
                return UnknownValue()
            if isinstance(rec, TaskResult):
                if field in rec.signature.outputs:
                    return UnknownValue()
            issues.append(f'Cannot resolve field access: {expr}')
            return UnknownValue()

        if expr.type == wf_node.NodeType.update:
            left = self._eval_expr(expr.left, imports, env, demanded, issues)
            right = self._eval_expr(expr.right, imports, env, demanded, issues)
            return self._merge_records(left, right)

        if expr.type == wf_node.NodeType.apply:
            fun_name = self._task_name(expr.fun)
            arg = self._eval_expr(expr.arg, imports, env, demanded, issues)
            if fun_name in imports:
                self._demand_task_inputs(fun_name, arg, imports, demanded)
                return TaskResult(fun_name, imports[fun_name].signature)
            return UnknownValue()

        if expr.type == wf_node.NodeType.block:
            local_env = dict(env)
            result = UnknownValue()
            for item in expr.body:
                if item.type == wf_node.NodeType.bind:
                    value = self._eval_expr(item.value, imports, local_env, demanded, issues)
                    local_env[item.id.name] = value
                    result = value
                else:
                    result = self._eval_expr(item, imports, local_env, demanded, issues)
            return result

        if expr.type == wf_node.NodeType.fun:
            return self._eval_function_body(expr, imports, env, demanded, issues)

        if expr.type == wf_node.NodeType.chain:
            left_name = self._task_name(expr.left)
            right_name = self._task_name(expr.right)
            if left_name in imports and right_name in imports:
                return TaskResult(right_name, imports[right_name].signature)
            return UnknownValue()

        return UnknownValue()

    def _merge_records(self, left, right):
        if isinstance(left, OpenRecord) and isinstance(right, OpenRecord):
            return OpenRecord(left.fields.union(right.fields))
        if isinstance(left, OpenRecord) and isinstance(right, ClosedRecord):
            return OpenRecord(left.fields.union(right.fields))
        if isinstance(left, ClosedRecord) and isinstance(right, OpenRecord):
            return OpenRecord(left.fields.union(right.fields))
        if isinstance(left, ClosedRecord) and isinstance(right, ClosedRecord):
            return ClosedRecord(left.fields.union(right.fields))
        if isinstance(left, TaskResult) and isinstance(right, ClosedRecord):
            return ClosedRecord(set(left.signature.outputs.keys()).union(right.fields))
        if isinstance(left, ClosedRecord) and isinstance(right, TaskResult):
            return ClosedRecord(left.fields.union(set(right.signature.outputs.keys())))
        if isinstance(left, OpenRecord) and isinstance(right, TaskResult):
            return OpenRecord(left.fields.union(set(right.signature.outputs.keys())))
        if isinstance(left, TaskResult) and isinstance(right, OpenRecord):
            return OpenRecord(set(left.signature.outputs.keys()).union(right.fields))
        if isinstance(left, TaskResult) and isinstance(right, TaskResult):
            return ClosedRecord(
                set(left.signature.outputs.keys()).union(
                    set(right.signature.outputs.keys())
                )
            )
        return UnknownValue()

    def _demand_task_inputs(self, task_name, arg, imports, demanded):
        signature = imports[task_name].signature
        available = self._available_fields(arg)
        for name, param in signature.inputs.items():
            if name in available:
                continue
            if param.type is not None and str(param.type.value).endswith('?'):
                continue
            demanded.add(name)

    def _available_fields(self, value):
        if isinstance(value, OpenRecord):
            return set(value.fields)
        if isinstance(value, ClosedRecord):
            return set(value.fields)
        if isinstance(value, TaskResult):
            return set(value.signature.outputs.keys())
        return set()

    def _task_name(self, expr):
        if expr.type == wf_node.NodeType.id:
            return expr.name
        if expr.type == wf_node.NodeType.chain:
            return self._task_name(expr.right)
        return None

    def _children(self, expr):
        if expr.type == wf_node.NodeType.block:
            return expr.body
        if expr.type == wf_node.NodeType.bind:
            return [expr.id, expr.value]
        if expr.type == wf_node.NodeType.rec:
            return list(expr.value.values())
        if expr.type == wf_node.NodeType.get:
            return [expr.rec, expr.member]
        if expr.type == wf_node.NodeType.update:
            return [expr.left, expr.right]
        if expr.type == wf_node.NodeType.fun:
            return [expr.param, expr.body]
        if expr.type == wf_node.NodeType.apply:
            return [expr.fun, expr.arg]
        if expr.type == wf_node.NodeType.chain:
            return [expr.left, expr.right]
        return []
