import os
from dataclasses import dataclass, field

from swl.semantic.task.type import Param, TaskSignature
from swl.semantic.wf import type as wf_type
from swl.semantic.wf.imports import load_import
from swl.semantic.wf.scope import chain_names
from swl.syntax.wf import builtins, node as wf_node


def _unknown_fields(fields=None):
    if isinstance(fields, dict):
        return dict(fields)
    return {name: UnknownValue() for name in (fields or [])}


@dataclass
class OpenRecord:
    fields: dict = field(default_factory=dict)

    def __init__(self, fields=None):
        self.fields = _unknown_fields(fields)


@dataclass
class ClosedRecord:
    fields: dict = field(default_factory=dict)

    def __init__(self, fields=None):
        self.fields = _unknown_fields(fields)


@dataclass
class FunctionValue:
    name: str
    signature: TaskSignature
    kind: str
    first_input: str = None
    param: str = None
    body: object = None
    env: dict = field(default_factory=dict)
    imports: dict = field(default_factory=dict)
    imported_check: object = None
    batch: bool = False

    def __post_init__(self):
        self.first_input = self.first_input or next(iter(self.signature.inputs.keys()), None)
        self.env = dict(self.env or {})
        self.imports = dict(self.imports or {})


@dataclass
class ClosureValue:
    function: FunctionValue
    bound_value: object = None

    def __post_init__(self):
        if self.bound_value is None:
            self.bound_value = ClosedRecord({})

    @property
    def signature(self):
        bound_fields = {
            name
            for name, value in getattr(self.bound_value, 'fields', {}).items()
            if is_explicitly_bound_value(value)
        }
        remaining = {
            name: param
            for name, param in self.function.signature.inputs.items()
            if name not in bound_fields
        }
        return TaskSignature(remaining, dict(self.function.signature.outputs), {})


@dataclass
class ComputationValue:
    function: FunctionValue
    arg_value: object = None
    output_value: object = None

    def __post_init__(self):
        if self.output_value is None:
            self.output_value = ClosedRecord({name: UnknownValue() for name in self.function.signature.outputs.keys()})

    @property
    def signature(self):
        return self.function.signature


@dataclass
class TableValue:
    columns: dict = field(default_factory=dict)
    placeholder: bool = False

    def __init__(self, columns=None, placeholder=False):
        self.columns = dict(columns or {})
        self.placeholder = placeholder


class UnknownValue:
    pass


@dataclass
class TypedValue:
    type: object = None


def is_explicitly_bound_value(value):
    return not isinstance(value, (UnknownValue, TypedValue))


def infer_inputs(checker, tree, imports):
    if tree.type != wf_node.NodeType.block or not tree.body:
        return set(), []

    final = tree.body[-1]
    issues = []
    demanded = set()

    if final.type == wf_node.NodeType.fun:
        env = eval_prefix_bindings(checker, tree.body[:-1], imports, demanded, issues)
        if uses_map(checker, final.body):
            table_demanded, table_issues, _ = eval_lambda_in_mode(checker, final, imports, env, 'table')
            if looks_like_stale_record_errors(table_issues):
                return table_demanded, []
            return table_demanded, table_issues

        env[final.param.name] = OpenRecord()
        eval_function_body(checker, final, imports, env, demanded, issues)
        return demanded, issues

    env = {}
    for expr in tree.body[:-1]:
        if expr.type == wf_node.NodeType.bind:
            if builtins.match_import(expr.value) is not None:
                continue
            env[expr.id.name] = eval_expr(checker, expr.value, imports, env, demanded, issues)
    value = eval_expr(checker, final, imports, env, demanded, issues)
    if not is_function_value(value):
        issues.append('Workflow must evaluate to a function')
    return demanded, issues


def _eval_block_items(checker, items, imports, env, demanded, issues):
    local_env = dict(env)
    result = UnknownValue()
    for expr in items:
        if expr.type == wf_node.NodeType.bind:
            value = eval_expr(checker, expr.value, imports, local_env, demanded, issues)
            local_env[expr.id.name] = value
            result = value
        else:
            result = eval_expr(checker, expr, imports, local_env, demanded, issues)
    return result


def eval_function_body(checker, fn, imports, env, demanded, issues):
    if fn.body.type != wf_node.NodeType.block:
        return eval_expr(checker, fn.body, imports, dict(env), demanded, issues)
    return _eval_block_items(checker, fn.body.body, imports, env, demanded, issues)


def _function_value_from_import(name, imported):
    return FunctionValue(
        name,
        imported.signature,
        imported.kind,
        imports=imported.check.imports if imported.check is not None else None,
        imported_check=imported.check,
        batch=imported.check.is_batch if imported.check is not None else False,
    )


def _inline_import_value(checker, path):
    base_dir = os.path.dirname(checker._loading[-1]) if checker._loading else '.'
    full_path = os.path.abspath(os.path.join(base_dir, path))
    stem = os.path.splitext(os.path.basename(full_path))[0]
    return _function_value_from_import(stem, load_import(checker, stem, full_path))


def _resolve_field_access(checker, rec, field, expr, demanded, issues):
    if isinstance(rec, OpenRecord):
        demanded.add(field)
        rec.fields[field] = UnknownValue()
        return rec.fields[field]
    if isinstance(rec, ClosedRecord):
        return rec.fields.get(field, UnknownValue())
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
        outputs = field_map(checker, rec.output_value)
        if outputs is not None and field in outputs:
            return outputs[field]
        if field in rec.signature.outputs:
            return UnknownValue()
    if isinstance(rec, (FunctionValue, ClosureValue)):
        issues.append(f'Cannot access field on function value: {expr}')
        return UnknownValue()
    issues.append(f'Cannot resolve field access: {expr}')
    return UnknownValue()


def eval_expr(checker, expr, imports, env, demanded, issues):
    if expr.type == wf_node.NodeType.id:
        if expr.name in env:
            return env[expr.name]
        if expr.name in imports:
            return _function_value_from_import(expr.name, imports[expr.name])
        return UnknownValue()

    if expr.type == wf_node.NodeType.num:
        return expr.value

    if expr.type == wf_node.NodeType.str:
        return expr.value

    if expr.type == wf_node.NodeType.rec:
        fields = {
            name: eval_expr(checker, value_expr, imports, env, demanded, issues)
            for name, value_expr in expr.value.items()
        }
        return ClosedRecord(fields)

    if expr.type == wf_node.NodeType.get:
        rec = eval_expr(checker, expr.rec, imports, env, demanded, issues)
        return _resolve_field_access(checker, rec, expr.member.name, expr, demanded, issues)

    if expr.type == wf_node.NodeType.update:
        left = eval_expr(checker, expr.left, imports, env, demanded, issues)
        right = eval_expr(checker, expr.right, imports, env, demanded, issues)
        return merge_update_values(checker, left, right, issues)

    if expr.type == wf_node.NodeType.apply:
        path = builtins.match_import(expr)
        if path is not None:
            return _inline_import_value(checker, path)
        map_parts = builtins.match_map(expr)
        map_by_parts = builtins.match_map_by(expr)
        if map_by_parts is not None:
            mapped_key = eval_expr(checker, map_by_parts[0], imports, env, demanded, issues)
            mapped_fun = eval_expr(checker, map_by_parts[1], imports, env, demanded, issues)
            mapped_arg = eval_expr(checker, map_by_parts[2], imports, env, demanded, issues)
            return apply_map_by(checker, mapped_fun, mapped_key, mapped_arg, issues)
        if map_parts is not None:
            mapped_fun = eval_expr(checker, map_parts[0], imports, env, demanded, issues)
            mapped_arg = eval_expr(checker, map_parts[1], imports, env, demanded, issues)
            return apply_map(checker, mapped_fun, mapped_arg, issues)
        if expr.fun.type == wf_node.NodeType.id and expr.fun.name == 'map':
            mapped_fun = eval_expr(checker, expr.arg, imports, env, demanded, issues)
            return apply_map_partial(checker, mapped_fun, issues)
        if expr.fun.type == wf_node.NodeType.id and expr.fun.name == 'map_by':
            mapped_fun = eval_expr(checker, expr.arg, imports, env, demanded, issues)
            return apply_map_by_partial(checker, mapped_fun, issues)
        fun = eval_expr(checker, expr.fun, imports, env, demanded, issues)
        arg = eval_expr(checker, expr.arg, imports, env, demanded, issues)
        if isinstance(fun, FunctionValue) and isinstance(arg, FunctionValue):
            if hasattr(checker, '_type_checker'):
                left_name = arg.name
                right_name = fun.name
                if left_name in checker._type_checker.signatures and right_name in checker._type_checker.signatures:
                    issues.extend(checker._type_checker.check_chain(left_name, right_name))
        result = apply(checker, fun, arg, demanded, issues)
        record_apply_satisfied(checker, expr, fun, arg, result)
        return result

    if expr.type == wf_node.NodeType.block:
        return _eval_block_items(checker, expr.body, imports, env, demanded, issues)

    if expr.type == wf_node.NodeType.fun:
        fn_issues = []
        fn_demanded = set()
        local_env = dict(env)
        local_env[expr.param.name] = OpenRecord()
        body_value = eval_function_body(checker, expr, imports, local_env, fn_demanded, fn_issues)
        batch = False
        if uses_map(checker, expr.body):
            _, table_issues, table_value = eval_lambda_in_mode(checker, expr, imports, env, 'table', fn_demanded)
            if value_is_table(table_value):
                body_value = table_value
                batch = True
                if isinstance(table_value, TableValue) and table_value.columns:
                    fn_demanded = set(table_value.columns.keys())
                else:
                    fn_demanded = {expr.param.name}
                fn_issues = [] if looks_like_stale_record_errors(table_issues) else table_issues
        body_outputs = signature_outputs(checker, body_value)

        if isinstance(body_value, (FunctionValue, ClosureValue)):
            body_outputs = dict(body_value.signature.outputs)
        if batch:
            root_input_type = infer_root_table_type(checker, expr, imports, env, fn_demanded)
            inputs = signature_inputs_from_table_type(checker, root_input_type)
            if not inputs:
                inputs = {expr.param.name: Param(expr.param.name, None)}
        else:
            if isinstance(body_value, ClosureValue):
                inputs = dict(body_value.signature.inputs)
            else:
                root_input_type = infer_root_record_type(checker, expr, imports, env, fn_demanded)
                inputs = signature_inputs_from_root_type(checker, root_input_type)
        signature = TaskSignature(inputs, body_outputs, {})
        return FunctionValue(
            expr.param.name,
            signature,
            'lambda',
            param=expr.param.name,
            body=expr.body,
            env=env,
            imports=imports,
            batch=batch or value_is_table(body_value),
        )

    if expr.type == wf_node.NodeType.chain:
        return eval_chain(checker, expr, imports)

    return UnknownValue()


def uses_map(checker, expr):
    if expr.type == wf_node.NodeType.apply and (builtins.match_map(expr) is not None or builtins.match_map_by(expr) is not None):
        return True
    from swl.semantic.wf.scope import children
    return any(uses_map(checker, child) for child in children(checker, expr))


def initial_param_value(mode, inferred_inputs=None):
    if mode == 'table':
        columns = {}
        if inferred_inputs:
            columns = {name: UnknownValue() for name in inferred_inputs}
        return TableValue(columns, placeholder=not inferred_inputs)
    if inferred_inputs is not None:
        return OpenRecord(inferred_inputs)
    return OpenRecord()


def eval_lambda_in_mode(checker, fn, imports, env, mode, inferred_inputs=None):
    local_env = dict(env)
    demanded = set()
    issues = []
    local_env[fn.param.name] = initial_param_value(mode, inferred_inputs)
    result = eval_function_body(checker, fn, imports, local_env, demanded, issues)
    if mode == 'table':
        table_value = local_env.get(fn.param.name)
        if isinstance(table_value, TableValue) and table_value.columns:
            demanded = set(table_value.columns.keys())
        else:
            demanded = {fn.param.name}
    return demanded, issues, result


def eval_prefix_bindings(checker, exprs, imports, demanded=None, issues=None):
    env = {}
    local_demanded = demanded if demanded is not None else set()
    local_issues = issues if issues is not None else []
    for expr in exprs:
        if expr.type != wf_node.NodeType.bind:
            continue
        if builtins.match_import(expr.value) is not None:
            continue
        env[expr.id.name] = eval_expr(checker, expr.value, imports, env, local_demanded, local_issues)
    return env


def looks_like_stale_record_errors(issues):
    return (
        ('map requires a tab argument' in issues or 'map_by requires a tab argument' in issues)
        and any(issue.startswith('Cannot resolve field access: (. ') for issue in issues)
    )


def value_is_table(value):
    return isinstance(value, TableValue)


def map_target(fun, issues):
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


def apply_map_partial(checker, fun, issues):
    target = map_target(fun, issues)
    if target is None:
        return UnknownValue()
    return FunctionValue(
        'map',
        TaskSignature({'f': Param('f', None), 'xs': Param('xs', None)}, dict(target.signature.outputs), {}),
        'builtin',
        first_input='f',
    )


def apply_map_by_partial(checker, fun, issues):
    target = map_target(fun, issues)
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


def _populate_table_inputs(arg, target):
    if not arg.columns:
        arg.columns.update({name: TypedValue(param.type) for name, param in target.signature.inputs.items()})
        arg.placeholder = False
        return
    for name, param in target.signature.inputs.items():
        if name not in arg.columns or isinstance(arg.columns[name], UnknownValue):
            arg.columns[name] = TypedValue(param.type)


def _apply_map_like(checker, fun, arg, issues, *, mode, key=None):
    target = map_target(fun, issues)
    if target is None:
        return UnknownValue()
    if mode == 'map_by' and not isinstance(key, str):
        issues.append('map_by requires a string key')
        return UnknownValue()
    if not isinstance(arg, TableValue):
        issues.append(f'{mode} requires a tab argument')
        return UnknownValue()
    _populate_table_inputs(arg, target)
    if mode == 'map_by':
        if key not in arg.columns:
            issues.append(f'Missing field on tab: {key}')
        if key not in target.signature.outputs:
            issues.append(f'map_by output must preserve grouping key: {key}')
    return TableValue({name: UnknownValue() for name in target.signature.outputs.keys()})


def apply_map(checker, fun, arg, issues):
    return _apply_map_like(checker, fun, arg, issues, mode='map')


def apply_map_by(checker, fun, key, arg, issues):
    return _apply_map_like(checker, fun, arg, issues, mode='map_by', key=key)


def record_apply_satisfied(checker, expr, fun, arg, result):
    if isinstance(result, ComputationValue):
        satisfied = set(result.function.signature.inputs.keys())
        checker._apply_satisfied[id(expr)] = satisfied
    elif isinstance(result, ClosureValue):
        bound = result.bound_value
        if isinstance(bound, (OpenRecord, ClosedRecord)):
            satisfied = {name for name, value in bound.fields.items() if is_explicitly_bound_value(value)}
        else:
            satisfied = set()
        checker._apply_satisfied[id(expr)] = satisfied


def _map_by_closure_result(checker, closure, arg, issues):
    bound_fields = dict(getattr(closure.bound_value, 'fields', {}))
    if 'f' in bound_fields and 'key' not in bound_fields:
        return ClosureValue(closure.function, ClosedRecord({'f': bound_fields['f'], 'key': arg}))
    if 'f' in bound_fields and 'key' in bound_fields:
        return apply_map_by(checker, bound_fields['f'], bound_fields['key'], arg, issues)
    return None


def _bound_application_value(checker, fun, arg):
    if isinstance(fun, ClosureValue):
        return merge_arg_values(checker, fun.bound_value, argument_value(checker, fun, arg))
    return argument_value(checker, fun, arg)


def apply(checker, fun, arg, demanded, issues):
    if not isinstance(fun, (ClosureValue, FunctionValue)):
        issues.append(f'Cannot apply non-function value: {fun}')
        return UnknownValue()
    if isinstance(fun, ClosureValue):
        return apply_closure(checker, fun, arg, demanded, issues)
    return apply_function(checker, fun, arg, demanded, issues)


def apply_closure(checker, closure, arg, demanded, issues=None):
    issues = issues or []
    if isinstance(closure.function, FunctionValue) and closure.function.kind == 'builtin' and closure.function.name == 'map_by':
        result = _map_by_closure_result(checker, closure, arg, issues)
        if result is not None:
            return result
    return application_result(checker, closure, _bound_application_value(checker, closure, arg), demanded, issues)


def apply_function(checker, fun, arg, demanded, issues=None):
    issues = issues or []
    if fun.kind == 'builtin' and fun.name == 'map_by':
        return apply_map_by_partial(checker, arg, issues)
    return application_result(checker, fun, _bound_application_value(checker, fun, arg), demanded, issues)


def application_result(checker, fun, bound, demanded, issues):
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
        return computation_value(checker, fun, bound)

    missing = required_missing(checker, signature, bound)
    for name in missing:
        demanded.add(name)
    if missing:
        return ClosureValue(base_fun, bound)
    if base_fun.kind == 'lambda' and base_fun.param is not None and base_fun.body is not None:
        local_env = dict(base_fun.env)
        local_env[base_fun.param] = bound
        return eval_expr(checker, base_fun.body, base_fun.imports, local_env, demanded, issues)
    if base_fun.kind == 'workflow' and base_fun.imported_check is not None:
        return apply_imported_workflow(checker, base_fun, bound, demanded, issues)
    return computation_value(checker, fun, bound)


def apply_imported_workflow(checker, fun, bound, demanded, issues):
    tree = fun.imported_check.tree
    if tree.type != wf_node.NodeType.block or not tree.body:
        return computation_value(checker, fun, bound)
    final = tree.body[-1]
    if final.type != wf_node.NodeType.fun:
        return computation_value(checker, fun, bound)
    env = {}
    for expr in tree.body[:-1]:
        if expr.type == wf_node.NodeType.bind:
            if builtins.match_import(expr.value) is not None:
                continue
            env[expr.id.name] = eval_expr(checker, expr.value, fun.imports, env, demanded, issues)
    env[final.param.name] = bound
    return eval_function_body(checker, final, fun.imports, env, demanded, issues)


def computation_value(checker, fun, bound):
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


def required_missing(checker, signature, available):
    missing = set()
    fields = field_map(checker, available)
    if fields is None:
        return set(signature.inputs.keys())
    for name, param in signature.inputs.items():
        if name in fields:
            continue
        if param.type is not None and str(param.type.value).endswith('?'):
            continue
        missing.add(name)
    return missing


def argument_value(checker, fun, arg):
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


def eval_chain(checker, expr, imports):
    names = chain_names(checker, expr)
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


def is_function_value(value):
    return isinstance(value, (FunctionValue, ClosureValue))


def merge_update_values(checker, left, right, issues):
    return merge_arg_values(checker, left, right)


def merge_tables(checker, left, right, issues):
    merged = dict(left.columns)
    for name, value in right.columns.items():
        if name in merged:
            left_type = value_type(checker, merged[name])
            right_type = value_type(checker, value)
            if left_type != wf_type.UNKNOWN and right_type != wf_type.UNKNOWN and left_type != right_type:
                issues.append(f'Type mismatch for "{name}": {left_type} -> {right_type}')
        merged[name] = value
    return TableValue(merged, placeholder=left.placeholder and right.placeholder)


def merge_arg_values(checker, left, right):
    left_fields = field_map(checker, left)
    right_fields = field_map(checker, right)
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
                left_type = value_type(checker, merged[name])
                right_type = value_type(checker, value)
                if left_type != wf_type.UNKNOWN and right_type != wf_type.UNKNOWN and left_type != right_type:
                    pass
                merged[name] = value
        else:
            merged[name] = value
    if isinstance(left, OpenRecord) or isinstance(right, OpenRecord):
        return OpenRecord(merged)
    return ClosedRecord(merged)


def field_map(checker, value):
    if isinstance(value, OpenRecord):
        return dict(value.fields)
    if isinstance(value, ClosedRecord):
        return dict(value.fields)
    if isinstance(value, FunctionValue):
        return {name: UnknownValue() for name in value.signature.outputs.keys()}
    if isinstance(value, ClosureValue):
        return {name: UnknownValue() for name in value.signature.outputs.keys()}
    if isinstance(value, ComputationValue):
        return field_map(checker, value.output_value)
    return None


def value_type(checker, value):
    if isinstance(value, TypedValue):
        return wf_type.scalar_from_name(getattr(value.type, 'value', None))
    if isinstance(value, TableValue):
        return wf_type.TableType({name: value_type(checker, item) for name, item in sorted(value.columns.items())})
    if isinstance(value, (OpenRecord, ClosedRecord)):
        return wf_type.RecordType({name: value_type(checker, item) for name, item in sorted(value.fields.items())}, open=isinstance(value, OpenRecord))
    if isinstance(value, ComputationValue):
        from swl.semantic.wf.signature import wf_function_type_from_signature
        return wf_function_type_from_signature(checker, value.signature, False).output
    if isinstance(value, (FunctionValue, ClosureValue)):
        from swl.semantic.wf.signature import wf_function_type_from_signature
        return wf_function_type_from_signature(checker, value.signature, getattr(getattr(value, 'function', value), 'batch', False))
    return wf_type.UNKNOWN


def signature_outputs(checker, value):
    if isinstance(value, ComputationValue):
        return dict(value.signature.outputs)
    if isinstance(value, ClosureValue):
        return dict(value.signature.outputs)
    if isinstance(value, FunctionValue):
        return dict(value.signature.outputs)
    if isinstance(value, ClosedRecord):
        return {
            name: param_from_value(checker, name, item)
            for name, item in sorted(value.fields.items())
        }
    if isinstance(value, OpenRecord):
        return {
            name: param_from_value(checker, name, item)
            for name, item in sorted(value.fields.items())
        }
    if isinstance(value, TableValue):
        return {
            name: param_from_value(checker, name, item)
            for name, item in sorted(value.columns.items())
        }
    return {}


def param_from_value(checker, name, value):
    typ = task_type_from_value(checker, value)
    return Param(name, typ)


def task_type_from_value(checker, value):
    if isinstance(value, TypedValue):
        return value.type
    if isinstance(value, ComputationValue) and len(value.signature.outputs) == 1:
        only = next(iter(value.signature.outputs.values()))
        return only.type
    if isinstance(value, (FunctionValue, ClosureValue, ComputationValue)):
        return None
    return None


def _typed_items(checker, items):
    return {name: value_type(checker, item) for name, item in sorted(items.items())}


def infer_root_record_type(checker, fn, imports, env, inferred_inputs):
    value = infer_root_input_value(checker, fn, imports, env, inferred_inputs, mode='record')
    fields = getattr(value, 'fields', {}) if isinstance(value, (OpenRecord, ClosedRecord)) else {}
    return wf_type.RecordType(_typed_items(checker, fields))


def infer_root_table_type(checker, fn, imports, env, inferred_inputs):
    value = infer_root_input_value(checker, fn, imports, env, inferred_inputs, mode='table')
    columns = getattr(value, 'columns', {}) if isinstance(value, TableValue) else {}
    if fn.param.name in columns and len(columns) > 1 and value_type(checker, columns[fn.param.name]) == wf_type.UNKNOWN:
        columns = {name: item for name, item in columns.items() if name != fn.param.name}
    return wf_type.TableType(_typed_items(checker, columns))


def infer_root_input_value(checker, fn, imports, env, inferred_inputs, mode):
    local_env = dict(env)
    local_env[fn.param.name] = initial_param_value(mode, inferred_inputs)
    eval_function_body(checker, fn, imports, local_env, set(), [])
    return local_env[fn.param.name]


def signature_inputs_from_root_type(checker, root_input_type):
    if isinstance(root_input_type, wf_type.RecordType):
        return {
            name: Param(name, task_type_from_wf_type(checker, typ))
            for name, typ in sorted(root_input_type.fields.items())
        }
    return {}


def signature_inputs_from_table_type(checker, root_input_type):
    if isinstance(root_input_type, wf_type.TableType):
        return {
            name: Param(name, task_type_from_wf_type(checker, typ))
            for name, typ in sorted(root_input_type.columns.items())
        }
    return {}


def task_type_from_wf_type(checker, typ):
    from swl.semantic.task.type import parse_type
    if isinstance(typ, wf_type.ScalarType) and typ.name in ('file', 'str', 'int', 'float'):
        return parse_type(typ.name)
    if isinstance(typ, wf_type.ArrayType) and isinstance(typ.item, wf_type.ScalarType):
        name = f'[{typ.item.name}]'
        return parse_type(name)
    return None


def inner_partial_remaining_inputs(checker, fn, imports, env):
    from swl.semantic.wf.scope import children
    body = getattr(fn, 'body', None)
    if getattr(body, 'type', None) != wf_node.NodeType.block or not getattr(body, 'body', None):
        return None
    result = body.body[-1]
    if getattr(result, 'type', None) != wf_node.NodeType.apply:
        return None
    if getattr(result.fun, 'type', None) != wf_node.NodeType.id:
        return None
    local_env = dict(env)
    local_env[fn.param.name] = initial_param_value('record', None)
    for expr in body.body[:-1]:
        if expr.type != wf_node.NodeType.bind:
            continue
        if builtins.match_import(expr.value) is not None:
            continue
        local_env[expr.id.name] = eval_expr(checker, expr.value, imports, local_env, set(), [])
    callee = local_env.get(result.fun.name)
    if not isinstance(callee, ClosureValue):
        return None
    if getattr(result.arg, 'type', None) != wf_node.NodeType.id or result.arg.name != fn.param.name:
        return None
    return callee.signature.inputs
