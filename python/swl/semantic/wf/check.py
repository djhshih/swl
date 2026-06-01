import os

from swl.semantic.task.type import Param, TaskSignature, TypeChecker, signature_from_task
from swl.semantic.wf import type as wf_type
from swl.syntax.task import bash as task_bash
from swl.syntax.task import interpolation as interp
from swl.syntax.task.parser import Parser as TaskParser
from swl.syntax.wf import builtins, node as wf_node
from swl.syntax.wf.parser import Parser as WfParser


class Import:
    def __init__(self, name: str, path: str, signature: TaskSignature, kind: str, check=None, task=None, parsed_body=None):
        self.name = name
        self.path = path
        self.signature = signature
        self.kind = kind
        self.check = check
        self.task = task
        self.parsed_body = parsed_body


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
    def __init__(self, name: str, signature: TaskSignature, kind: str, first_input: str = None, param=None, body=None, env=None, imports=None, imported_check=None, batch=False):
        self.name = name
        self.signature = signature
        self.kind = kind
        self.first_input = first_input or self._first_input_name()
        self.param = param
        self.body = body
        self.env = dict(env or {})
        self.imports = dict(imports or {})
        self.imported_check = imported_check
        self.batch = batch

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
        bound_fields = set()
        if isinstance(self.bound_value, (OpenRecord, ClosedRecord)):
            for name, value in self.bound_value.fields.items():
                if _is_explicitly_bound_value(value):
                    bound_fields.add(name)
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


class TableValue:
    def __init__(self, columns=None, placeholder=False):
        self.columns = dict(columns or {})
        self.placeholder = placeholder


class UnknownValue:
    pass


class TypedValue:
    def __init__(self, typ=None):
        self.type = typ


def _is_explicitly_bound_value(value):
    return not isinstance(value, (UnknownValue, TypedValue))


def _bash_interp_vars(word):
    vars = []
    parts = word.parts if isinstance(word, interp.Word) else [word]
    for part in parts:
        if isinstance(part, interp.Var):
            vars.append(part.name)
    return vars


def _validate_bash_variables(parsed_body, input_names, context_name):
    errors = []
    defined = set(input_names)
    for stmt in parsed_body.statements:
        if isinstance(stmt, task_bash.Assignment):
            refs = _bash_interp_vars(stmt.value)
            for ref in refs:
                if ref not in defined:
                    errors.append(
                        f'Unresolved variable "${ref}" in assignment to "{stmt.name}" '
                        f'in {context_name}'
                    )
            defined.add(stmt.name)
        elif isinstance(stmt, task_bash.Command):
            for word in stmt.words:
                for ref in _bash_interp_vars(word):
                    if ref not in defined:
                        errors.append(
                            f'Unresolved variable "${ref}" in command "{stmt.text[:60]}" '
                            f'in {context_name}'
                        )
    return errors


class WorkflowCheck:
    def __init__(self, tree, imports, errors, inferred_inputs, signature=None, is_batch=False, workflow_type=None, root_input_type=None, root_output_type=None):
        self.tree = tree
        self.imports = imports
        self.errors = list(errors)
        self.inferred_inputs = inferred_inputs
        self.signature = signature
        self.is_batch = is_batch
        self.workflow_type = workflow_type
        self.root_input_type = root_input_type
        self.root_output_type = root_output_type

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
        self._apply_satisfied = {}

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
            self._type_checker = checker
            errors.extend(self._check_chains(tree, checker))
            inferred_inputs, infer_errors = self._infer_inputs(tree, imports)
            signature, is_batch, workflow_type, root_input_type, root_output_type = self._build_workflow_signature(tree, imports, inferred_inputs, errors)
            errors.extend(infer_errors)
            if signature is None:
                errors.append('Workflow must evaluate to a function')
            return WorkflowCheck(tree, imports, errors, inferred_inputs, signature, is_batch, workflow_type, root_input_type, root_output_type)
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
            path = builtins.match_import(expr.value)
            if path is None:
                continue
            full_path = os.path.join(base_dir, path)
            imports[expr.id.name] = self._load_import(expr.id.name, full_path)
        return imports

    def _load_import(self, name: str, path: str) -> Import:
        if path.endswith('.sh'):
            src = self._read_file(path)
            task = TaskParser().parse(src)
            parsed_body = task_bash.Parser().parse(task.body)
            signature = signature_from_task(task)
            input_names = set(signature.inputs.keys())
            var_errors = _validate_bash_variables(parsed_body, input_names, f'task "{name}" ({path})')
            if var_errors:
                raise ValueError('\n'.join(var_errors))
            return Import(name, path, signature, 'task', task=task, parsed_body=parsed_body)
        if path.endswith('.swl'):
            check = self.load(path)
            if check.signature is None:
                raise ValueError(f'Imported workflow does not produce a signature: {path}')
            return Import(name, path, check.signature, 'workflow', check=check)
        raise ValueError(f'Unrecognized import path: {path}')

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
            env = self._eval_prefix_bindings(tree.body[:-1], imports, demanded, issues)
            if self._uses_map(final.body):
                table_demanded, table_issues, _ = self._eval_lambda_in_mode(final, imports, env, 'table')
                if self._looks_like_stale_record_errors(table_issues):
                    return table_demanded, []
                return table_demanded, table_issues

            env[final.param.name] = OpenRecord()
            self._eval_function_body(final, imports, env, demanded, issues)
            return demanded, issues

        env = {}
        for expr in tree.body[:-1]:
            if expr.type == wf_node.NodeType.bind:
                if builtins.match_import(expr.value) is not None:
                    continue
                env[expr.id.name] = self._eval_expr(expr.value, imports, env, demanded, issues)
        value = self._eval_expr(final, imports, env, demanded, issues)
        if not self._is_function_value(value):
            issues.append('Workflow must evaluate to a function')
        return demanded, issues

    def _build_workflow_signature(self, tree, imports, inferred_inputs, issues):
        if tree.type != wf_node.NodeType.block or not tree.body:
            return None, False, None, None, None

        final = tree.body[-1]
        if final.type == wf_node.NodeType.fun:
            env = self._eval_prefix_bindings(tree.body[:-1], imports, set(), issues)

            if self._uses_map(final.body):
                demanded, mode_issues, chosen_result = self._eval_lambda_in_mode(final, imports, env, 'table', inferred_inputs)
                if self._value_is_table(chosen_result):
                    issues[:] = [
                        issue for issue in issues
                        if not issue.startswith('Cannot resolve field access: (. ')
                        and issue not in ('map requires a tab argument', 'map_by requires a tab argument')
                    ]
                    root_input_type = self._infer_root_table_type(final, imports, env, inferred_inputs)
                    root_output_type = self._wf_output_type(chosen_result, table_context=True)
                    inputs = self._signature_inputs_from_table_type(root_input_type)
                    if not inputs:
                        inputs = {final.param.name: Param(final.param.name, None)}
                    outputs = self._signature_outputs(chosen_result)
                    signature = TaskSignature(inputs, outputs, {})
                    return signature, True, wf_type.FunctionType(root_input_type, root_output_type), root_input_type, root_output_type
                filtered_issues = [] if self._looks_like_stale_record_errors(mode_issues) else mode_issues
                issues.extend(filtered_issues)
                root_input_type = self._infer_root_table_type(final, imports, env, inferred_inputs)
                root_output_type = self._wf_output_type(chosen_result, table_context=True)
                inputs = self._signature_inputs_from_table_type(root_input_type)
                if not inputs:
                    inputs = {final.param.name: Param(final.param.name, None)}
                outputs = self._signature_outputs(chosen_result)
                signature = TaskSignature(inputs, outputs, {})
                return signature, True, wf_type.FunctionType(root_input_type, root_output_type), root_input_type, root_output_type

            record_env = dict(env)
            record_env[final.param.name] = self._initial_param_value('record', inferred_inputs)
            chosen_result = self._eval_function_body(final, imports, record_env, set(), list(issues))
            root_output_type = self._wf_output_type(chosen_result, table_context=False)
            partial_inputs = self._inner_partial_remaining_inputs(final, imports, env)
            if isinstance(chosen_result, ClosureValue):
                inputs = dict(chosen_result.signature.inputs)
                root_input_type = wf_type.RecordType({
                    name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
                    for name, param in sorted(inputs.items())
                }, open=True)
            elif partial_inputs is not None:
                inputs = dict(partial_inputs)
                root_input_type = wf_type.RecordType({
                    name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
                    for name, param in sorted(inputs.items())
                }, open=True)
            else:
                root_input_type = self._infer_root_record_type(final, imports, env, inferred_inputs)
                inputs = self._signature_inputs_from_root_type(root_input_type)
            outputs = self._signature_outputs(chosen_result)
            signature = TaskSignature(inputs, outputs, {})
            return signature, False, wf_type.FunctionType(root_input_type, root_output_type), root_input_type, root_output_type

        env = {}
        for expr in tree.body[:-1]:
            if expr.type == wf_node.NodeType.bind:
                if builtins.match_import(expr.value) is not None:
                    continue
                env[expr.id.name] = self._eval_expr(expr.value, imports, env, set(), issues)
        result = self._eval_expr(final, imports, env, set(), issues)
        if isinstance(result, FunctionValue):
            wft = self._wf_function_type_from_signature(result.signature, result.batch)
            return result.signature, result.batch, wft, getattr(wft, 'input', None), getattr(wft, 'output', None)
        if isinstance(result, ClosureValue):
            signature = result.signature
            wft = self._wf_function_type_from_signature(signature, result.function.batch)
            return signature, result.function.batch, wft, getattr(wft, 'input', None), getattr(wft, 'output', None)
        return None, False, None, None, None

    def _signature_outputs(self, value):
        if isinstance(value, ComputationValue):
            return dict(value.signature.outputs)
        if isinstance(value, ClosureValue):
            return dict(value.signature.outputs)
        if isinstance(value, FunctionValue):
            return dict(value.signature.outputs)
        if isinstance(value, ClosedRecord):
            return {
                name: self._param_from_value(name, item)
                for name, item in sorted(value.fields.items())
            }
        if isinstance(value, OpenRecord):
            return {
                name: self._param_from_value(name, item)
                for name, item in sorted(value.fields.items())
            }
        if isinstance(value, TableValue):
            return {
                name: self._param_from_value(name, item)
                for name, item in sorted(value.columns.items())
            }
        return {}

    def _wf_output_type(self, value, table_context=False):
        if isinstance(value, TableValue):
            return wf_type.TableType({name: wf_type.UNKNOWN for name in sorted(value.columns.keys())})
        if isinstance(value, (ClosedRecord, OpenRecord)):
            fields = getattr(value, 'fields', {})
            return wf_type.RecordType({name: wf_type.UNKNOWN for name in sorted(fields.keys())}, open=isinstance(value, OpenRecord))
        if isinstance(value, (ComputationValue, ClosureValue, FunctionValue)):
            outputs = self._signature_outputs(value)
            return wf_type.RecordType({name: wf_type.UNKNOWN for name in sorted(outputs.keys())}, open=False)
        return wf_type.RecordType({}, open=True)

    def _wf_function_type_from_signature(self, signature, is_batch):
        if signature is None:
            return None
        output = wf_type.RecordType({
            name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
            for name, param in sorted(signature.outputs.items())
        }, open=False)
        if is_batch:
            input_type = wf_type.TableType({
                name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
                for name, param in sorted(signature.inputs.items())
            })
        else:
            input_type = wf_type.RecordType({
                name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
                for name, param in sorted(signature.inputs.items())
            }, open=True)
        return wf_type.FunctionType(input_type, output)

    def _signature_inputs_from_root_type(self, root_input_type):
        if isinstance(root_input_type, wf_type.RecordType):
            return {
                name: Param(name, self._task_type_from_wf_type(typ))
                for name, typ in sorted(root_input_type.fields.items())
            }
        return {}

    def _signature_inputs_from_table_type(self, root_input_type):
        if isinstance(root_input_type, wf_type.TableType):
            return {
                name: Param(name, self._task_type_from_wf_type(typ))
                for name, typ in sorted(root_input_type.columns.items())
            }
        return {}

    def _task_type_from_wf_type(self, typ):
        from swl.semantic.task.type import parse_type
        if isinstance(typ, wf_type.ScalarType) and typ.name in ('file', 'str', 'int', 'float'):
            return parse_type(typ.name)
        if isinstance(typ, wf_type.ArrayType) and isinstance(typ.item, wf_type.ScalarType):
            name = f'[{typ.item.name}]'
            return parse_type(name)
        return None

    def _infer_root_table_type(self, fn, imports, env, inferred_inputs):
        value = self._infer_root_input_value(fn, imports, env, inferred_inputs, mode='table')
        columns = getattr(value, 'columns', {}) if isinstance(value, TableValue) else {}
        if fn.param.name in columns and len(columns) > 1 and self._value_type(columns[fn.param.name]) == wf_type.UNKNOWN:
            columns = {name: item for name, item in columns.items() if name != fn.param.name}
        return wf_type.TableType({name: self._value_type(item) for name, item in sorted(columns.items())})

    def _infer_root_input_value(self, fn, imports, env, inferred_inputs, mode):
        local_env = dict(env)
        local_env[fn.param.name] = self._initial_param_value(mode, inferred_inputs)
        self._eval_function_body(fn, imports, local_env, set(), [])
        return local_env[fn.param.name]

    def _inner_partial_remaining_inputs(self, fn, imports, env):
        body = getattr(fn, 'body', None)
        if getattr(body, 'type', None) != wf_node.NodeType.block or not getattr(body, 'body', None):
            return None
        result = body.body[-1]
        if getattr(result, 'type', None) != wf_node.NodeType.apply:
            return None
        if getattr(result.fun, 'type', None) != wf_node.NodeType.id:
            return None
        local_env = dict(env)
        local_env[fn.param.name] = self._initial_param_value('record', None)
        for expr in body.body[:-1]:
            if expr.type != wf_node.NodeType.bind:
                continue
            if builtins.match_import(expr.value) is not None:
                continue
            local_env[expr.id.name] = self._eval_expr(expr.value, imports, local_env, set(), [])
        callee = local_env.get(result.fun.name)
        if not isinstance(callee, ClosureValue):
            return None
        if getattr(result.arg, 'type', None) != wf_node.NodeType.id or result.arg.name != fn.param.name:
            return None
        return callee.signature.inputs

    def _value_type(self, value):
        if isinstance(value, TypedValue):
            return wf_type.scalar_from_name(getattr(value.type, 'value', None))
        if isinstance(value, TableValue):
            return wf_type.TableType({name: self._value_type(item) for name, item in sorted(value.columns.items())})
        if isinstance(value, (OpenRecord, ClosedRecord)):
            return wf_type.RecordType({name: self._value_type(item) for name, item in sorted(value.fields.items())}, open=isinstance(value, OpenRecord))
        if isinstance(value, ComputationValue):
            return self._wf_function_type_from_signature(value.signature, False).output
        if isinstance(value, (FunctionValue, ClosureValue)):
            return self._wf_function_type_from_signature(value.signature, getattr(getattr(value, 'function', value), 'batch', False))
        return wf_type.UNKNOWN

    def _infer_root_record_type(self, fn, imports, env, inferred_inputs):
        local_env = dict(env)
        local_env[fn.param.name] = self._initial_param_value('record', inferred_inputs)
        result = self._eval_function_body(fn, imports, local_env, set(), [])
        if isinstance(result, ClosureValue):
            return wf_type.RecordType({
                name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
                for name, param in sorted(result.signature.inputs.items())
            }, open=True)
        if isinstance(result, ComputationValue) and isinstance(result.arg_value, OpenRecord):
            fields = result.arg_value.fields
            return wf_type.RecordType({name: self._value_type(item) for name, item in sorted(fields.items())}, open=True)
        value = local_env[fn.param.name]
        fields = getattr(value, 'fields', {}) if isinstance(value, (OpenRecord, ClosedRecord)) else {}
        return wf_type.RecordType({name: self._value_type(item) for name, item in sorted(fields.items())}, open=True)

    def _task_type_from_value(self, value):
        if isinstance(value, TypedValue):
            return value.type
        if isinstance(value, ComputationValue) and len(value.signature.outputs) == 1:
            only = next(iter(value.signature.outputs.values()))
            return only.type
        if isinstance(value, (FunctionValue, ClosureValue, ComputationValue)):
            return None
        return None

    def _param_from_value(self, name, value):
        typ = self._task_type_from_value(value)
        return Param(name, typ)

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
                    batch=imported.check.is_batch if imported.check is not None else False,
                )
            return UnknownValue()

        if expr.type == wf_node.NodeType.num:
            return expr.value

        if expr.type == wf_node.NodeType.str:
            return expr.value

        if expr.type == wf_node.NodeType.rec:
            fields = {
                name: self._eval_expr(value_expr, imports, env, demanded, issues)
                for name, value_expr in expr.value.items()
            }
            return ClosedRecord(fields)

        if expr.type == wf_node.NodeType.get:
            rec_expr = expr.rec
            rec = self._eval_expr(rec_expr, imports, env, demanded, issues)
            field = expr.member.name
            if isinstance(rec, OpenRecord):
                demanded.add(field)
                rec.fields[field] = UnknownValue()
                return rec.fields[field]
            if isinstance(rec, ClosedRecord):
                if field in rec.fields:
                    return rec.fields[field]
                return UnknownValue()
            if isinstance(rec, TableValue):
                if field in rec.columns:
                    return rec.columns[field]
                if rec.placeholder:
                    issues.append('map requires a tab argument')
                    return UnknownValue()
                issues.append(f'Missing field on tab: {field}')
                rec.columns[field] = UnknownValue()
                return rec.columns[field]
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
            return self._merge_update_values(left, right, issues)

        if expr.type == wf_node.NodeType.apply:
            path = builtins.match_import(expr)
            if path is not None:
                base_dir = os.path.dirname(self._loading[-1]) if self._loading else '.'
                full_path = os.path.abspath(os.path.join(base_dir, path))
                stem = os.path.splitext(os.path.basename(full_path))[0]
                imp = self._load_import(stem, full_path)
                return FunctionValue(
                    stem,
                    imp.signature,
                    imp.kind,
                    imports=imp.check.imports if imp.check is not None else None,
                    imported_check=imp.check,
                    batch=imp.check.is_batch if imp.check is not None else False,
                )
            map_parts = builtins.match_map(expr)
            map_by_parts = builtins.match_map_by(expr)
            if map_by_parts is not None:
                mapped_key = self._eval_expr(map_by_parts[0], imports, env, demanded, issues)
                mapped_fun = self._eval_expr(map_by_parts[1], imports, env, demanded, issues)
                mapped_arg = self._eval_expr(map_by_parts[2], imports, env, demanded, issues)
                return self._apply_map_by(mapped_fun, mapped_key, mapped_arg, issues)
            if map_parts is not None:
                mapped_fun = self._eval_expr(map_parts[0], imports, env, demanded, issues)
                mapped_arg = self._eval_expr(map_parts[1], imports, env, demanded, issues)
                return self._apply_map(mapped_fun, mapped_arg, issues)
            if expr.fun.type == wf_node.NodeType.id and expr.fun.name == 'map':
                mapped_fun = self._eval_expr(expr.arg, imports, env, demanded, issues)
                return self._apply_map_partial(mapped_fun, issues)
            if expr.fun.type == wf_node.NodeType.id and expr.fun.name == 'map_by':
                mapped_fun = self._eval_expr(expr.arg, imports, env, demanded, issues)
                return self._apply_map_by_partial(mapped_fun, issues)
            fun = self._eval_expr(expr.fun, imports, env, demanded, issues)
            arg = self._eval_expr(expr.arg, imports, env, demanded, issues)
            if isinstance(fun, FunctionValue) and isinstance(arg, FunctionValue):
                if hasattr(self, '_type_checker'):
                    left_name = arg.name
                    right_name = fun.name
                    if left_name in self._type_checker.signatures and right_name in self._type_checker.signatures:
                        issues.extend(self._type_checker.check_chain(left_name, right_name))
            result = self._apply(fun, arg, demanded, issues)
            self._record_apply_satisfied(expr, fun, arg, result)
            return result

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
            batch = False
            if self._uses_map(expr.body):
                _, table_issues, table_value = self._eval_lambda_in_mode(expr, imports, env, 'table', fn_demanded)
                if self._value_is_table(table_value):
                    body_value = table_value
                    batch = True
                    if isinstance(table_value, TableValue) and table_value.columns:
                        fn_demanded = set(table_value.columns.keys())
                    else:
                        fn_demanded = {expr.param.name}
                    fn_issues = [] if self._looks_like_stale_record_errors(table_issues) else table_issues
            body_outputs = self._signature_outputs(body_value)

            if isinstance(body_value, (FunctionValue, ClosureValue)):
                body_outputs = dict(body_value.signature.outputs)
            if batch:
                root_input_type = self._infer_root_table_type(expr, imports, env, fn_demanded)
                inputs = self._signature_inputs_from_table_type(root_input_type)
                if not inputs:
                    inputs = {expr.param.name: Param(expr.param.name, None)}
            else:
                if isinstance(body_value, ClosureValue):
                    inputs = dict(body_value.signature.inputs)
                else:
                    root_input_type = self._infer_root_record_type(expr, imports, env, fn_demanded)
                    inputs = self._signature_inputs_from_root_type(root_input_type)
            signature = TaskSignature(inputs, body_outputs, {})
            return FunctionValue(
                expr.param.name,
                signature,
                'lambda',
                param=expr.param.name,
                body=expr.body,
                env=env,
                imports=imports,
                batch=batch or self._value_is_table(body_value),
            )

        if expr.type == wf_node.NodeType.chain:
            return self._eval_chain(expr, imports)

        return UnknownValue()

    def _uses_map(self, expr):
        if expr.type == wf_node.NodeType.apply and (builtins.match_map(expr) is not None or builtins.match_map_by(expr) is not None):
            return True
        return any(self._uses_map(child) for child in self._children(expr))

    def _initial_param_value(self, mode, inferred_inputs=None):
        if mode == 'table':
            columns = {}
            if inferred_inputs:
                columns = {name: UnknownValue() for name in inferred_inputs}
            return TableValue(columns, placeholder=not inferred_inputs)
        if inferred_inputs is not None:
            return OpenRecord(inferred_inputs)
        return OpenRecord()

    def _eval_lambda_in_mode(self, fn, imports, env, mode, inferred_inputs=None):
        local_env = dict(env)
        demanded = set()
        issues = []
        local_env[fn.param.name] = self._initial_param_value(mode, inferred_inputs)
        result = self._eval_function_body(fn, imports, local_env, demanded, issues)
        if mode == 'table':
            table_value = local_env.get(fn.param.name)
            if isinstance(table_value, TableValue) and table_value.columns:
                demanded = set(table_value.columns.keys())
            else:
                demanded = {fn.param.name}
        return demanded, issues, result

    def _eval_prefix_bindings(self, exprs, imports, demanded=None, issues=None):
        env = {}
        local_demanded = demanded if demanded is not None else set()
        local_issues = issues if issues is not None else []
        for expr in exprs:
            if expr.type != wf_node.NodeType.bind:
                continue
            if builtins.match_import(expr.value) is not None:
                continue
            env[expr.id.name] = self._eval_expr(expr.value, imports, env, local_demanded, local_issues)
        return env

    def _looks_like_stale_record_errors(self, issues):
        return (
            ('map requires a tab argument' in issues or 'map_by requires a tab argument' in issues)
            and any(issue.startswith('Cannot resolve field access: (. ') for issue in issues)
        )

    def _value_is_table(self, value):
        return isinstance(value, TableValue)

    def _map_target(self, fun, issues):
        target = fun.function if isinstance(fun, ClosureValue) else fun
        if not isinstance(target, FunctionValue):
            issues.append('map requires a function value')
            return None
        if target.batch:
            issues.append('map on batch workflow is not supported')
            return None
        if len(target.signature.outputs) == 0:
            issues.append('map requires a function with rec output')
            return None
        return target

    def _apply_map_partial(self, fun, issues):
        target = self._map_target(fun, issues)
        if target is None:
            return UnknownValue()
        return FunctionValue(
            'map',
            TaskSignature({'f': Param('f', None), 'xs': Param('xs', None)}, dict(target.signature.outputs), {}),
            'builtin',
            first_input='f',
        )

    def _apply_map_by_partial(self, fun, issues):
        target = self._map_target(fun, issues)
        if target is None:
            return UnknownValue()
        return ClosureValue(
            FunctionValue(
                'map_by',
                TaskSignature({'key': Param('key', None), 'xs': Param('xs', None)}, dict(target.signature.outputs), {}),
                'builtin',
                first_input='key',
            ),
            ClosedRecord({'f': target}),
        )

    def _apply_map(self, fun, arg, issues):
        target = self._map_target(fun, issues)
        if target is None:
            return UnknownValue()
        if not isinstance(arg, TableValue):
            issues.append('map requires a tab argument')
            return UnknownValue()
        if not arg.columns:
            arg.columns.update({
                name: TypedValue(param.type)
                for name, param in target.signature.inputs.items()
            })
            arg.placeholder = False
        else:
            for name, param in target.signature.inputs.items():
                if name not in arg.columns or isinstance(arg.columns[name], UnknownValue):
                    arg.columns[name] = TypedValue(param.type)
        return TableValue({name: UnknownValue() for name in target.signature.outputs.keys()})

    def _apply_map_by(self, fun, key, arg, issues):
        target = self._map_target(fun, issues)
        if target is None:
            return UnknownValue()
        if not isinstance(key, str):
            issues.append('map_by requires a string key')
            return UnknownValue()
        if not isinstance(arg, TableValue):
            issues.append('map_by requires a tab argument')
            return UnknownValue()
        if not arg.columns:
            arg.columns.update({
                name: TypedValue(param.type)
                for name, param in target.signature.inputs.items()
            })
            arg.placeholder = False
        else:
            for name, param in target.signature.inputs.items():
                if name not in arg.columns or isinstance(arg.columns[name], UnknownValue):
                    arg.columns[name] = TypedValue(param.type)
        if key not in arg.columns:
            issues.append(f'Missing field on tab: {key}')
        if key not in target.signature.outputs:
            issues.append(f'map_by output must preserve grouping key: {key}')
        return TableValue({name: UnknownValue() for name in target.signature.outputs.keys()})

    def _record_apply_satisfied(self, expr, fun, arg, result):
        if isinstance(result, ComputationValue):
            satisfied = set(result.function.signature.inputs.keys())
            self._apply_satisfied[id(expr)] = satisfied
        elif isinstance(result, ClosureValue):
            bound = result.bound_value
            if isinstance(bound, (OpenRecord, ClosedRecord)):
                satisfied = {name for name, value in bound.fields.items() if _is_explicitly_bound_value(value)}
            else:
                satisfied = set()
            self._apply_satisfied[id(expr)] = satisfied

    def _apply(self, fun, arg, demanded, issues):
        if isinstance(fun, ClosureValue):
            return self._apply_closure(fun, arg, demanded, issues)
        if isinstance(fun, FunctionValue):
            return self._apply_function(fun, arg, demanded, issues)
        issues.append(f'Cannot apply non-function value: {fun}')
        return UnknownValue()

    def _apply_closure(self, closure, arg, demanded, issues=None):
        if isinstance(closure.function, FunctionValue) and closure.function.kind == 'builtin' and closure.function.name == 'map_by':
            bound_fields = dict(getattr(closure.bound_value, 'fields', {}))
            if 'f' in bound_fields and 'key' not in bound_fields:
                return ClosureValue(closure.function, ClosedRecord({'f': bound_fields['f'], 'key': arg}))
            if 'f' in bound_fields and 'key' in bound_fields:
                return self._apply_map_by(bound_fields['f'], bound_fields['key'], arg, issues or [])
        bound = self._merge_arg_values(closure.bound_value, self._argument_value(closure, arg))
        return self._application_result(closure, bound, demanded, issues or [])

    def _apply_function(self, fun, arg, demanded, issues=None):
        if isinstance(fun, FunctionValue) and fun.kind == 'builtin' and fun.name == 'map_by':
            return self._apply_map_by_partial(arg, issues or [])
        bound = self._argument_value(fun, arg)
        return self._application_result(fun, bound, demanded, issues or [])

    def _application_result(self, fun, bound, demanded, issues):
        signature = fun.signature
        base_fun = fun.function if isinstance(fun, ClosureValue) else fun
        if isinstance(bound, OpenRecord):
            for name, param in signature.inputs.items():
                if name in bound.fields:
                    current = bound.fields[name]
                    if isinstance(current, (UnknownValue, TypedValue)):
                        bound.fields[name] = TypedValue(param.type)
                    continue
                if param.type is not None and str(param.type.value).endswith('?'):
                    continue
                demanded.add(name)
                bound.fields[name] = TypedValue(param.type)
            return self._computation_value(fun, bound)

        missing = self._required_missing(signature, bound)
        for name in missing:
            demanded.add(name)
        if missing:
            return ClosureValue(base_fun, bound)
        if base_fun.kind == 'lambda' and base_fun.param is not None and base_fun.body is not None:
            local_env = dict(base_fun.env)
            local_env[base_fun.param] = bound
            return self._eval_expr(base_fun.body, base_fun.imports, local_env, demanded, issues)
        if base_fun.kind == 'workflow' and base_fun.imported_check is not None:
            return self._apply_imported_workflow(base_fun, bound, demanded, issues)
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
                if builtins.match_import(expr.value) is not None:
                    continue
                env[expr.id.name] = self._eval_expr(expr.value, fun.imports, env, demanded, issues)
        env[final.param.name] = bound
        return self._eval_function_body(final, fun.imports, env, demanded, issues)

    def _computation_value(self, fun, bound):
        effective_fun = fun
        if isinstance(fun, ClosureValue):
            effective_fun = FunctionValue(
                fun.function.name,
                fun.signature,
                fun.function.kind,
                param=fun.function.param,
                body=fun.function.body,
                env=fun.function.env,
                imports=fun.function.imports,
                imported_check=fun.function.imported_check,
                batch=fun.function.batch,
            )
        return ComputationValue(
            effective_fun,
            bound,
            ClosedRecord({name: UnknownValue() for name in effective_fun.signature.outputs.keys()}),
        )

    def _required_missing(self, signature, available):
        missing = set()
        fields = self._field_map(available)
        if fields is None:
            return set(signature.inputs.keys())
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
            if arg.kind == 'lambda':
                return arg
            return ClosedRecord({name: UnknownValue() for name in arg.signature.outputs.keys()})
        if isinstance(arg, TableValue):
            return arg
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

    def _merge_update_values(self, left, right, issues):
        if isinstance(left, TableValue) and isinstance(right, TableValue):
            return self._merge_tables(left, right, issues)
        if isinstance(left, TableValue):
            return self._merge_table_record(left, right)
        if isinstance(right, TableValue):
            return self._merge_table_record(right, left)
        return self._merge_records(left, right)

    def _merge_table_record(self, table, record):
        record_fields = self._field_map(record)
        if record_fields is None:
            return UnknownValue()
        merged = dict(table.columns)
        merged.update(record_fields)
        return TableValue(merged, placeholder=table.placeholder)

    def _merge_tables(self, left, right, issues):
        merged = dict(left.columns)
        for name, value in right.columns.items():
            if name in merged:
                left_type = self._value_type(merged[name])
                right_type = self._value_type(value)
                if left_type != wf_type.UNKNOWN and right_type != wf_type.UNKNOWN and left_type != right_type:
                    issues.append(f'Type mismatch for "{name}": {left_type} -> {right_type}')
            merged[name] = value
        return TableValue(merged, placeholder=left.placeholder and right.placeholder)

    def _merge_arg_values(self, left, right):
        left_fields = self._field_map(left)
        right_fields = self._field_map(right)
        if left_fields is None or right_fields is None:
            return UnknownValue()
        merged = dict(left_fields)
        for name, value in right_fields.items():
            if name in merged:
                if isinstance(merged[name], UnknownValue) and isinstance(value, TypedValue):
                    merged[name] = value
                elif isinstance(value, UnknownValue) and isinstance(merged[name], TypedValue):
                    pass
                else:
                    left_type = self._value_type(merged[name])
                    right_type = self._value_type(value)
                    if left_type != wf_type.UNKNOWN and right_type != wf_type.UNKNOWN and left_type != right_type:
                        pass
                    merged[name] = value
            else:
                merged[name] = value
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
