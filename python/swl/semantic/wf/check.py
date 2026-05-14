import os

from swl.semantic.task.type import Param, TaskSignature, TypeChecker, signature_from_task
from swl.syntax.task.parser import Parser as TaskParser
from swl.syntax.wf import node as wf_node
from swl.syntax.wf.parser import Parser as WfParser


class Import:
    def __init__(self, name: str, path: str, signature: TaskSignature, kind: str, check=None):
        self.name = name
        self.path = path
        self.signature = signature
        self.kind = kind
        self.check = check


class OpenRecord:
    def __init__(self, fields=None):
        if isinstance(fields, dict):
            self.fields = dict(fields)
        else:
            self.fields = {name: UnknownValue() for name in (fields or [])}


class ClosedRecord:
    def __init__(self, fields=None):
        if isinstance(fields, dict):
            self.fields = dict(fields)
        else:
            self.fields = {name: UnknownValue() for name in (fields or [])}


class FunctionValue:
    def __init__(self, name: str, signature: TaskSignature, kind: str, first_input: str = None, param=None, body=None, env=None, imports=None, imported_check=None):
        self.name = name
        self.signature = signature
        self.kind = kind
        self.first_input = first_input or self._first_input_name()
        self.param = param
        self.body = body
        self.env = dict(env or {})
        self.imports = dict(imports or {})
        self.imported_check = imported_check

    def _first_input_name(self):
        for name in self.signature.inputs.keys():
            return name
        return None


class ClosureValue:
    def __init__(self, function: FunctionValue, bound_value=None):
        self.function = function
        self.bound_value = bound_value or ClosedRecord({})

    @property
    def signature(self):
        remaining = {}
        bound_fields = set(self.bound_value.fields.keys()) if isinstance(self.bound_value, (OpenRecord, ClosedRecord)) else set()
        for name, param in self.function.signature.inputs.items():
            if name not in bound_fields:
                remaining[name] = param
        return TaskSignature(remaining, dict(self.function.signature.outputs), {})


class ComputationValue:
    def __init__(self, function: FunctionValue, arg_value=None, output_value=None):
        self.function = function
        self.arg_value = arg_value
        self.output_value = output_value or ClosedRecord({
            name: UnknownValue() for name in function.signature.outputs.keys()
        })

    @property
    def signature(self):
        return self.function.signature


class UnknownValue:
    pass


class WorkflowCheck:
    def __init__(self, tree, imports, errors, inferred_inputs, signature=None):
        self.tree = tree
        self.imports = imports
        self.errors = list(errors)
        self.inferred_inputs = inferred_inputs
        self.signature = signature

    @property
    def chain_errors(self):
        return self.errors

    @property
    def issues(self):
        return self.errors


class Checker:
    def __init__(self, files=None):
        self._loading = []
        self.files = files or {}

    def load(self, path: str) -> WorkflowCheck:
        full_path = os.path.abspath(path)
        if full_path in self._loading:
            raise ValueError(f'Circular workflow import: {full_path}')
        self._loading.append(full_path)
        try:
            src = self._read_file(full_path)
            tree = WfParser().parse(src)
            imports = self._load_imports(tree, os.path.dirname(full_path))
            errors = self._check_scope(tree)
            checker = TypeChecker()
            for imported in imports.values():
                checker.add_task(imported.name, imported.signature)
            errors.extend(self._check_chains(tree, checker))
            inferred_inputs, infer_errors = self._infer_inputs(tree, imports)
            errors.extend(infer_errors)
            signature = self._build_workflow_signature(tree, imports, inferred_inputs, errors)
            if signature is None:
                errors.append('Workflow must evaluate to a function')
            return WorkflowCheck(tree, imports, errors, inferred_inputs, signature)
        finally:
            self._loading.pop()

    def load_content(self, content: str, path: str = '<memory>.swl') -> WorkflowCheck:
        full_path = os.path.abspath(path)
        old = self.files.get(full_path)
        self.files[full_path] = content
        try:
            return self.load(full_path)
        finally:
            if old is None:
                del self.files[full_path]
            else:
                self.files[full_path] = old

    def _read_file(self, path: str) -> str:
        if path in self.files:
            return self.files[path]
        with open(path, 'r') as f:
            return f.read()

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
            task = TaskParser().parse(self._read_file(path))
            return Import(name, path, signature_from_task(task), 'task')
        if path.endswith('.swl'):
            check = self.load(path)
            if check.signature is None:
                raise ValueError(f'Imported workflow does not produce a signature: {path}')
            return Import(name, path, check.signature, 'workflow', check=check)
        raise ValueError(f'Unrecognized import path: {path}')

    def _match_import(self, expr):
        if expr.type != wf_node.NodeType.apply:
            return None
        if expr.fun.type != wf_node.NodeType.id or expr.fun.name != 'import':
            return None
        if expr.arg.type != wf_node.NodeType.str:
            return None
        return expr.arg.value

    def _check_scope(self, tree):
        errors = []
        self._walk_scope(tree, set(), errors)
        return errors

    def _walk_scope(self, expr, scope, errors):
        if expr.type == wf_node.NodeType.block:
            local = set(scope)
            for item in expr.body:
                if item.type == wf_node.NodeType.bind:
                    if item.id.name in local:
                        errors.append(f'Duplicate binding in scope: {item.id.name}')
                    else:
                        local.add(item.id.name)
                    self._walk_scope(item.value, local, errors)
                else:
                    self._walk_scope(item, local, errors)
            return

        if expr.type == wf_node.NodeType.fun:
            local = {expr.param.name}
            self._walk_scope(expr.body, local, errors)
            return

        for child in self._children(expr):
            self._walk_scope(child, scope, errors)

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
            env = {}
            for expr in tree.body[:-1]:
                if expr.type == wf_node.NodeType.bind:
                    if self._match_import(expr.value) is not None:
                        continue
                    env[expr.id.name] = self._eval_expr(expr.value, imports, env, demanded, issues)
            env[final.param.name] = OpenRecord()
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
            env = {}
            for expr in tree.body[:-1]:
                if expr.type == wf_node.NodeType.bind:
                    if self._match_import(expr.value) is not None:
                        continue
                    env[expr.id.name] = self._eval_expr(expr.value, imports, env, set(), issues)
            env[final.param.name] = OpenRecord(inferred_inputs)
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
            return {name: Param(name, None) for name in sorted(value.fields.keys())}
        if isinstance(value, OpenRecord):
            return {name: Param(name, None) for name in sorted(value.fields.keys())}
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
                return FunctionValue(
                    expr.name,
                    imported.signature,
                    imported.kind,
                    imports=imported.check.imports if imported.check is not None else None,
                    imported_check=imported.check,
                )
            return UnknownValue()

        if expr.type == wf_node.NodeType.rec:
            fields = {
                name: self._eval_expr(value_expr, imports, env, demanded, issues)
                for name, value_expr in expr.value.items()
            }
            return ClosedRecord(fields)

        if expr.type == wf_node.NodeType.get:
            rec = self._eval_expr(expr.rec, imports, env, demanded, issues)
            field = expr.member.name
            if isinstance(rec, OpenRecord):
                demanded.add(field)
                rec.fields[field] = UnknownValue()
                return rec.fields[field]
            if isinstance(rec, ClosedRecord):
                if field in rec.fields:
                    return rec.fields[field]
                return UnknownValue()
            if isinstance(rec, ComputationValue):
                outputs = self._field_map(rec.output_value)
                if outputs is not None and field in outputs:
                    return outputs[field]
                if field in rec.signature.outputs:
                    return UnknownValue()
            if isinstance(rec, (FunctionValue, ClosureValue)):
                issues.append(f'Cannot access field on function value: {expr}')
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
            fn_issues = []
            fn_demanded = set()
            local_env = dict(env)
            local_env[expr.param.name] = OpenRecord()
            body_value = self._eval_function_body(expr, imports, local_env, fn_demanded, fn_issues)
            body_outputs = self._signature_outputs(body_value)
            inputs = {
                name: Param(name, None)
                for name in sorted(fn_demanded)
            }
            signature = TaskSignature(inputs, body_outputs, {})
            return FunctionValue(
                expr.param.name,
                signature,
                'lambda',
                param=expr.param.name,
                body=expr.body,
                env=env,
                imports=imports,
            )

        if expr.type == wf_node.NodeType.chain:
            return self._eval_chain(expr, imports)

        return UnknownValue()

    def _apply(self, fun, arg, demanded, issues):
        if isinstance(fun, ClosureValue):
            return self._apply_closure(fun, arg, demanded, issues)
        if isinstance(fun, FunctionValue):
            return self._apply_function(fun, arg, demanded, issues)
        issues.append(f'Cannot apply non-function value: {fun}')
        return UnknownValue()

    def _apply_closure(self, closure, arg, demanded, issues=None):
        bound = self._merge_arg_values(closure.bound_value, self._argument_value(closure.function, arg))
        return self._application_result(closure.function, bound, demanded, issues or [])

    def _apply_function(self, fun, arg, demanded, issues=None):
        bound = self._argument_value(fun, arg)
        return self._application_result(fun, bound, demanded, issues or [])

    def _application_result(self, fun, bound, demanded, issues):
        if isinstance(bound, OpenRecord):
            for name, param in fun.signature.inputs.items():
                if name in bound.fields:
                    continue
                if param.type is not None and str(param.type.value).endswith('?'):
                    continue
                demanded.add(name)
                bound.fields[name] = UnknownValue()
            return self._computation_value(fun, bound)

        missing = self._required_missing(fun.signature, bound)
        for name in missing:
            demanded.add(name)
        if missing:
            return ClosureValue(fun, bound)
        if fun.kind == 'lambda' and fun.param is not None and fun.body is not None:
            local_env = dict(fun.env)
            local_env[fun.param] = bound
            return self._eval_expr(fun.body, fun.imports, local_env, demanded, issues)
        if fun.kind == 'workflow' and fun.imported_check is not None:
            return self._apply_imported_workflow(fun, bound, demanded, issues)
        return self._computation_value(fun, bound)

    def _apply_imported_workflow(self, fun, bound, demanded, issues):
        tree = fun.imported_check.tree
        if tree.type != wf_node.NodeType.block or not tree.body:
            return self._computation_value(fun, bound)
        final = tree.body[-1]
        if final.type != wf_node.NodeType.fun:
            return self._computation_value(fun, bound)
        env = {}
        for expr in tree.body[:-1]:
            if expr.type == wf_node.NodeType.bind:
                if self._match_import(expr.value) is not None:
                    continue
                env[expr.id.name] = self._eval_expr(expr.value, fun.imports, env, demanded, issues)
        env[final.param.name] = bound
        return self._eval_function_body(final, fun.imports, env, demanded, issues)

    def _computation_value(self, fun, bound):
        return ComputationValue(
            fun,
            bound,
            ClosedRecord({name: UnknownValue() for name in fun.signature.outputs.keys()}),
        )

    def _required_missing(self, signature, available):
        missing = set()
        fields = self._field_map(available)
        for name, param in signature.inputs.items():
            if name in fields:
                continue
            if param.type is not None and str(param.type.value).endswith('?'):
                continue
            missing.add(name)
        return missing

    def _argument_value(self, fun, arg):
        if isinstance(arg, (OpenRecord, ClosedRecord)):
            return arg
        if isinstance(arg, ComputationValue):
            return ClosedRecord({name: UnknownValue() for name in arg.signature.outputs.keys()})
        if isinstance(arg, ClosureValue):
            return ClosedRecord({name: UnknownValue() for name in arg.signature.outputs.keys()})
        if isinstance(arg, FunctionValue):
            return ClosedRecord({name: UnknownValue() for name in arg.signature.outputs.keys()})
        if fun.first_input is not None:
            return ClosedRecord({fun.first_input: arg})
        return ClosedRecord({})

    def _eval_chain(self, expr, imports):
        names = self._chain_names(expr)
        if not names:
            return UnknownValue()
        missing = [name for name in names if name not in imports]
        if missing:
            return UnknownValue()

        signatures = [imports[name].signature for name in names]
        inputs = dict(signatures[0].inputs)
        available = set(signatures[0].inputs.keys()).union(signatures[0].outputs.keys())
        outputs = {}
        run = {}
        for sig in signatures:
            for name, param in sig.inputs.items():
                if name not in available:
                    inputs[name] = param
            available.update(sig.inputs.keys())
            available.update(sig.outputs.keys())
            outputs.update(sig.outputs)
            run.update(sig.run)
        kind = imports[names[-1]].kind
        return FunctionValue(names[-1], TaskSignature(inputs, outputs, run), kind)

    def _is_function_value(self, value):
        return isinstance(value, (FunctionValue, ClosureValue))

    def _merge_records(self, left, right):
        return self._merge_arg_values(left, right)

    def _merge_arg_values(self, left, right):
        left_fields = self._field_map(left)
        right_fields = self._field_map(right)
        if left_fields is None or right_fields is None:
            return UnknownValue()
        merged = dict(left_fields)
        merged.update(right_fields)
        if isinstance(left, OpenRecord) or isinstance(right, OpenRecord):
            return OpenRecord(merged)
        return ClosedRecord(merged)

    def _field_map(self, value):
        if isinstance(value, OpenRecord):
            return dict(value.fields)
        if isinstance(value, ClosedRecord):
            return dict(value.fields)
        if isinstance(value, FunctionValue):
            return {name: UnknownValue() for name in value.signature.outputs.keys()}
        if isinstance(value, ClosureValue):
            return {name: UnknownValue() for name in value.signature.outputs.keys()}
        if isinstance(value, ComputationValue):
            return self._field_map(value.output_value)
        return None

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
