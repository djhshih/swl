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


class FunctionValue:
    def __init__(self, name: str, signature: TaskSignature, kind: str, first_input: str = None):
        self.name = name
        self.signature = signature
        self.kind = kind
        self.first_input = first_input or self._first_input_name()

    def _first_input_name(self):
        for name in self.signature.inputs.keys():
            return name
        return None


class ClosureValue:
    def __init__(self, function: FunctionValue, bound_fields=None):
        self.function = function
        self.bound_fields = set(bound_fields or [])

    @property
    def signature(self):
        remaining = {}
        for name, param in self.function.signature.inputs.items():
            if name not in self.bound_fields:
                remaining[name] = param
        return TaskSignature(remaining, dict(self.function.signature.outputs), {})


class ComputationValue:
    def __init__(self, function: FunctionValue, available_fields=None):
        self.function = function
        self.available_fields = set(available_fields or [])

    @property
    def signature(self):
        return self.function.signature


class UnknownValue:
    pass


class WorkflowCheck:
    def __init__(self, tree, imports, errors, inferred_inputs, signature=None):
        self.tree = tree
        self.imports = imports
        self.chain_errors = errors
        self.errors = errors
        self.inferred_inputs = inferred_inputs
        self.issues = errors
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
            errors = self._check_chains(tree, checker)
            inferred_inputs, infer_errors = self._infer_inputs(tree, imports)
            errors.extend(infer_errors)
            signature = self._build_workflow_signature(tree, imports, inferred_inputs, errors)
            if signature is None:
                errors.append('Workflow must evaluate to a function')
            return WorkflowCheck(tree, imports, errors, inferred_inputs, signature)
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
                if self._match_import(expr.value) is not None:
                    continue
                env[expr.id.name] = self._eval_expr(expr.value, imports, env, demanded, issues)
        value = self._eval_expr(final, imports, env, demanded, issues)
        if not self._is_function_value(value):
            issues.append('Workflow must evaluate to a function')
        return demanded, issues

    def _build_workflow_signature(self, tree, imports, inferred_inputs, issues):
        if tree.type != wf_node.NodeType.block or not tree.body:
            return None

        final = tree.body[-1]
        if final.type == wf_node.NodeType.fun:
            env = {final.param.name: OpenRecord(inferred_inputs)}
            result = self._eval_function_body(final, imports, env, set(), issues)
            inputs = {
                name: Param(name, None)
                for name in sorted(inferred_inputs)
            }
            outputs = self._signature_outputs(result)
            return TaskSignature(inputs, outputs, {})

        env = {}
        for expr in tree.body[:-1]:
            if expr.type == wf_node.NodeType.bind:
                if self._match_import(expr.value) is not None:
                    continue
                env[expr.id.name] = self._eval_expr(expr.value, imports, env, set(), issues)
        result = self._eval_expr(final, imports, env, set(), issues)
        if isinstance(result, FunctionValue):
            return result.signature
        if isinstance(result, ClosureValue):
            return result.signature
        return None

    def _signature_outputs(self, value):
        if isinstance(value, ComputationValue):
            return dict(value.signature.outputs)
        if isinstance(value, ClosureValue):
            return dict(value.signature.outputs)
        if isinstance(value, FunctionValue):
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
            if expr.name in env:
                return env[expr.name]
            if expr.name in imports:
                imported = imports[expr.name]
                return FunctionValue(expr.name, imported.signature, imported.kind)
            return UnknownValue()

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
            if isinstance(rec, ComputationValue):
                if field in rec.signature.outputs:
                    return UnknownValue()
            if isinstance(rec, ClosureValue):
                if field in rec.signature.outputs:
                    return UnknownValue()
            issues.append(f'Cannot resolve field access: {expr}')
            return UnknownValue()

        if expr.type == wf_node.NodeType.update:
            left = self._eval_expr(expr.left, imports, env, demanded, issues)
            right = self._eval_expr(expr.right, imports, env, demanded, issues)
            return self._merge_records(left, right)

        if expr.type == wf_node.NodeType.apply:
            path = self._match_import(expr)
            if path is not None:
                return UnknownValue()
            fun = self._eval_expr(expr.fun, imports, env, demanded, issues)
            arg = self._eval_expr(expr.arg, imports, env, demanded, issues)
            return self._apply(fun, arg, demanded, issues)

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
            return UnknownValue()

        if expr.type == wf_node.NodeType.chain:
            return self._eval_chain(expr, imports)

        return UnknownValue()

    def _apply(self, fun, arg, demanded, issues):
        if isinstance(fun, ClosureValue):
            return self._apply_closure(fun, arg, demanded)
        if isinstance(fun, FunctionValue):
            return self._apply_function(fun, arg, demanded)
        issues.append(f'Cannot apply non-function value: {fun}')
        return UnknownValue()

    def _apply_closure(self, closure, arg, demanded):
        available = set(closure.bound_fields)
        available.update(self._available_fields_for_apply(closure.function, arg))
        return self._application_result(closure.function, available, demanded)

    def _apply_function(self, fun, arg, demanded):
        available = self._available_fields_for_apply(fun, arg)
        return self._application_result(fun, available, demanded)

    def _application_result(self, fun, available, demanded):
        if isinstance(available, OpenRecord):
            for name, param in fun.signature.inputs.items():
                if name in available.fields:
                    continue
                if param.type is not None and str(param.type.value).endswith('?'):
                    continue
                demanded.add(name)
                available.fields.add(name)
            return ComputationValue(fun, available.fields)

        missing = self._required_missing(fun.signature, available)
        for name in missing:
            demanded.add(name)
        if missing:
            return ClosureValue(fun, available)
        return ComputationValue(fun, available)

    def _required_missing(self, signature, available):
        missing = set()
        for name, param in signature.inputs.items():
            if name in available:
                continue
            if param.type is not None and str(param.type.value).endswith('?'):
                continue
            missing.add(name)
        return missing

    def _available_fields_for_apply(self, fun, arg):
        if isinstance(arg, OpenRecord):
            return arg
        if isinstance(arg, (ClosedRecord, ComputationValue, ClosureValue, FunctionValue)):
            return self._available_fields(arg)
        if fun.first_input is not None:
            return {fun.first_input}
        return set()

    def _eval_chain(self, expr, imports):
        names = self._chain_names(expr)
        if not names:
            return UnknownValue()
        missing = [name for name in names if name not in imports]
        if missing:
            return UnknownValue()

        first = imports[names[0]].signature
        outputs = dict(first.outputs)
        for name in names[1:]:
            outputs.update(imports[name].signature.outputs)
        kind = imports[names[-1]].kind
        return FunctionValue(names[-1], TaskSignature(first.inputs, outputs, {}), kind)

    def _is_function_value(self, value):
        return isinstance(value, (FunctionValue, ClosureValue))

    def _merge_records(self, left, right):
        if isinstance(left, OpenRecord) and isinstance(right, OpenRecord):
            return OpenRecord(left.fields.union(right.fields))
        if isinstance(left, OpenRecord) and isinstance(right, ClosedRecord):
            return OpenRecord(left.fields.union(right.fields))
        if isinstance(left, ClosedRecord) and isinstance(right, OpenRecord):
            return OpenRecord(left.fields.union(right.fields))
        if isinstance(left, ClosedRecord) and isinstance(right, ClosedRecord):
            return ClosedRecord(left.fields.union(right.fields))
        if isinstance(left, ComputationValue) and isinstance(right, ClosedRecord):
            return ClosedRecord(set(left.signature.outputs.keys()).union(right.fields))
        if isinstance(left, ClosedRecord) and isinstance(right, ComputationValue):
            return ClosedRecord(left.fields.union(set(right.signature.outputs.keys())))
        if isinstance(left, OpenRecord) and isinstance(right, ComputationValue):
            return OpenRecord(left.fields.union(set(right.signature.outputs.keys())))
        if isinstance(left, ComputationValue) and isinstance(right, OpenRecord):
            return OpenRecord(set(left.signature.outputs.keys()).union(right.fields))
        if isinstance(left, ComputationValue) and isinstance(right, ComputationValue):
            return ClosedRecord(
                set(left.signature.outputs.keys()).union(
                    set(right.signature.outputs.keys())
                )
            )
        return UnknownValue()

    def _available_fields(self, value):
        if isinstance(value, OpenRecord):
            return set(value.fields)
        if isinstance(value, ClosedRecord):
            return set(value.fields)
        if isinstance(value, FunctionValue):
            return set(value.signature.outputs.keys())
        if isinstance(value, ClosureValue):
            return set(value.signature.outputs.keys())
        if isinstance(value, ComputationValue):
            return set(value.signature.outputs.keys())
        return set()

    def _chain_names(self, expr):
        if expr.type == wf_node.NodeType.id:
            return [expr.name]
        if expr.type == wf_node.NodeType.chain:
            return self._chain_names(expr.left) + self._chain_names(expr.right)
        return []

    def _task_name(self, expr):
        if expr.type == wf_node.NodeType.id:
            return expr.name
        if expr.type == wf_node.NodeType.chain:
            names = self._chain_names(expr)
            if names:
                return names[-1]
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
