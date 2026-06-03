from swl.ir import node as ir
from swl.dag.context import ForceEnv
from swl.dag.evaluator import _available_inputs, _force_map, _value_key
from swl.dag.node import Field, ForcedFunction, Input, Literal, Record, StepCall
from swl.types import normalize_swl_type


def _input(forcer, name, spec=None):
    if name not in forcer.inputs:
        typ = None
        desc = None
        if spec is not None:
            typ = spec.type.value if getattr(spec, 'type', None) is not None else None
            desc = spec.desc
        optional = bool(typ and typ.endswith('?'))
        forcer.inputs[name] = Input(name, typ, desc, optional)
    return forcer.inputs[name]


def _function_signature(forcer, function):
    if isinstance(function, ir.Function):
        return function.signature
    if isinstance(function, ir.Lambda):
        return function.signature
    return None


def _forced_signature(forcer, value):
    signature = value.signature or _function_signature(forcer, value.function)
    if signature is None:
        return None
    if value.bound is None:
        return signature
    available = _available_inputs(forcer, value.bound)
    inputs = {
        name: param
        for name, param in signature.inputs.items()
        if name not in available
    }
    return type(signature)(inputs, dict(signature.outputs), dict(signature.run))


def _force_root(forcer, value):
    if not isinstance(value, ForcedFunction):
        return value
    signature = _forced_signature(forcer, value)
    if signature is None:
        return value
    if _is_unsaturated_builtin_map(forcer, value, signature):
        return _materialize_partial_map_workflow(forcer, value, signature)
    if _is_batch_function(forcer, value):
        table_name = getattr(getattr(value.function, 'body', None), 'param', None) or getattr(value.function, 'param', None) or 'xs'
        for name, spec in signature.inputs.items():
            _input(forcer, name, spec)
        from swl.dag.evaluator import _apply
        return _apply(forcer, value, Record({table_name: _input(forcer, table_name)}))
    arg = Record({name: _input(forcer, name, spec) for name, spec in signature.inputs.items()})
    from swl.dag.evaluator import _apply
    return _apply(forcer, value, arg)


def _is_batch_function(forcer, fn):
    function = fn.function if isinstance(fn, ForcedFunction) else fn
    return bool(getattr(function, 'is_batch', False))


def _is_unsaturated_builtin_map(forcer, value, signature):
    return (
        isinstance(value, ForcedFunction)
        and isinstance(value.function, ir.Function)
        and value.function.kind == 'builtin'
        and value.function.name in ('map', 'map_by')
        and signature is not None
        and 'xs' in signature.inputs
    )


def _materialize_partial_map_workflow(forcer, value, signature):
    if not isinstance(value.bound, Record) or 'f' not in value.bound.fields:
        raise ValueError('partial map workflow is missing bound function')
    fn = value.bound.fields['f']
    if not isinstance(fn, ForcedFunction):
        raise ValueError(f'partial map requires a normalized function, got: {fn!r}')
    xs_spec = signature.inputs['xs']
    xs = _input(forcer, 'xs', xs_spec)
    key = None
    if isinstance(value.function, ir.Function) and value.function.name == 'map_by':
        if 'key' not in value.bound.fields:
            raise ValueError('partial map_by workflow is missing bound grouping key')
        key_value = value.bound.fields['key']
        key = key_value.value if isinstance(key_value, Literal) else None
    mapped = _force_map(forcer, fn, xs, key=key)
    if isinstance(mapped, StepCall):
        return Record({name: Field(mapped, name) for name in mapped.outputs})
    return mapped


def _materialize_workflow_dag(forcer, function):
    from swl.dag.forcer import make_force_state
    sub_forcer = make_force_state(files=forcer.lowerer.checker.loader)
    return _materialize_workflow_dag_impl(sub_forcer, function)


def _materialize_workflow_dag_impl(forcer, function):
    signature = function.signature
    if signature is None:
        from swl.dag.forcer import force
        return force(function.body, state=forcer)
    arg = Record({name: _input(forcer, name, spec) for name, spec in signature.inputs.items()})
    if isinstance(function.body, ir.Lambda):
        from swl.dag.evaluator import _apply
        value = _apply(forcer, ForcedFunction(function.body, None, signature), arg)
    else:
        from swl.dag.evaluator import force_value
        value = force_value(forcer, function.body, ForceEnv())
        if isinstance(value, ForcedFunction):
            from swl.dag.evaluator import _apply
            value = _apply(forcer, value, arg)
    if isinstance(value, ForcedFunction):
        raise ValueError(f'Workflow materialization did not resolve to a concrete value: {function.path}')
    from swl.dag.finalize import _finalize_dag
    return _finalize_dag(forcer, value)


def _tool_definition(forcer, path):
    if path in forcer.tool_defs:
        return forcer.tool_defs[path]
    cached = forcer.lowerer.checker.loader.get_parsed_task(path)
    if cached is not None:
        task, signature, parsed_body = cached
    else:
        from swl.semantic.task.type import signature_from_task
        from swl.semantic.wf.bashvars import _validate_bash_variables
        from swl.semantic.wf.imports import read_file
        from swl.syntax.task import bash as task_bash
        from swl.syntax.task.parser import Parser as TaskParser
        src = read_file(forcer.lowerer.checker, path)
        task = TaskParser().parse(src)
        signature = signature_from_task(task)
        parsed_body = task_bash.Parser().parse(task.body)
        known_vars = set(signature.inputs.keys()) | set(signature.run.keys())
        var_errors = _validate_bash_variables(parsed_body, known_vars, f'task at "{path}"')
        if var_errors:
            raise ValueError('\n'.join(var_errors))
    definition = {
        'doc': task.annotation.doc,
        'body': task.body,
        'inputs': _build_task_section(forcer, task, signature, 'in'),
        'outputs': _build_task_section(forcer, task, signature, 'out'),
        'run': _build_task_section(forcer, task, signature, 'run'),
    }
    forcer.tool_defs[path] = definition
    return definition


def _workflow_definition(forcer, function):
    key = ('workflow', function.path)
    if key in forcer.tool_defs:
        return forcer.tool_defs[key]
    body_dag = _materialize_workflow_dag(forcer, function).to_dict()
    definition = {
        'class': 'Workflow',
        'dag': body_dag,
        'inputs': {name: {'type': normalize_swl_type(spec.type.value if getattr(spec, 'type', None) is not None else None), 'desc': spec.desc} for name, spec in function.signature.inputs.items()},
        'outputs': {name: {'type': normalize_swl_type(spec.type.value if getattr(spec, 'type', None) is not None else None)} for name, spec in function.signature.outputs.items()},
        'body': '',
        'run': {},
    }
    forcer.tool_defs[key] = definition
    return definition


def _build_task_section(forcer, task, signature, kind):
    built = {}
    for section in task.annotation.sections:
        if section.kind.value != kind:
            continue
        for param in section.params:
            for name in param.names:
                built[name] = _build_task_param(forcer, signature, kind, name, param)
    return built


def _build_task_param(forcer, signature, kind, name, param):
    if kind == 'run':
        semantic = signature.run[name]
        built = {
            'type': semantic.type.value if semantic.type is not None else None,
            'value': semantic.parsed_default,
        }
        if param.desc is not None:
            built['desc'] = param.desc
        return built
    semantic = signature.inputs.get(name) if kind == 'in' else signature.outputs.get(name)
    built = {
        'type': semantic.type.value if semantic and semantic.type is not None else param.type,
    }
    if param.desc is not None:
        built['desc'] = param.desc
    if param.default is not None:
        built['default'] = _interp_to_dict(param.default)
        if kind == 'out':
            _validate_output_default_glob(forcer, name, param.default)
    return built


def _validate_output_default_glob(forcer, name, default):
    text = _interp_word_to_text(default)
    if '[' in text or ']' in text:
        depth = 0
        for c in text:
            if c == '[':
                depth += 1
            elif c == ']':
                depth -= 1
            if depth < 0 or depth > 1:
                raise ValueError(f'Malformed glob pattern in output "{name}" default: {text!r}')
        if depth != 0:
            raise ValueError(f'Malformed glob pattern in output "{name}" default: {text!r}')
    if '*' in text or '?' in text:
        # Glob metacharacters are allowed; only bracket structure is validated here.
        return


def _interp_to_dict(value):
    from swl.syntax.task import interpolation as interp
    if isinstance(value, interp.Word):
        return {'kind': 'word', 'parts': [_interp_to_dict(part) for part in value.parts]}
    if isinstance(value, interp.Literal):
        return {'kind': 'literal', 'text': value.text}
    if isinstance(value, interp.Var):
        return {'kind': 'var', 'name': value.name}
    if isinstance(value, interp.Expr):
        return {'kind': 'expr', 'text': value.text}
    return None


def _interp_word_to_text(value):
    from swl.syntax.task import interpolation as interp
    if isinstance(value, interp.Word):
        return ''.join(_interp_word_to_text(part) for part in value.parts)
    if isinstance(value, interp.Literal):
        return value.text
    if isinstance(value, interp.Var):
        return '${' + value.name + '}'
    if isinstance(value, interp.Expr):
        return '${' + value.text + '}'
    return ''
