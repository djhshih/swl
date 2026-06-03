from swl.dag.binding import binding_to_dict
from swl.dag.context import _SENTINEL
from swl.dag.evaluator import _can_project_record_fields, _project_field, _project_input
from swl.dag.merge import _value_key
from swl.dag.node import Field, Input, Merge, Record, StepCall, TableSource
from swl.types import normalize_swl_type


def _step_output_record(call, outputs):
    return Record({name: Field(call, name) for name in outputs})


def _emit_step_call(forcer, function, bound, satisfied, *, run, task, step_type='task'):
    inputs = _normalize_task_inputs(forcer, function, bound)
    _validate_bindings(forcer, inputs, satisfied, function.name)
    key = _step_call_key(forcer, function, inputs)
    if key in forcer.step_cache:
        return forcer.step_cache[key]
    if step_type == 'task':
        forcer.call_counter += 1
    outputs = list(function.signature.outputs.keys())
    call = StepCall(
        id=_step_id(forcer, function.name),
        path=function.path,
        bindings=inputs,
        outputs=outputs,
        run=run,
        task=task,
        deps=sorted(_step_dependencies(forcer, inputs)),
        type=step_type,
    )
    forcer.steps.append(call)
    result = _step_output_record(call, outputs)
    forcer.step_cache[key] = result
    return result


def _emit_task_call(forcer, function, bound, satisfied=None):
    from swl.dag.tooldefs import _tool_definition
    return _emit_step_call(
        forcer,
        function,
        bound,
        satisfied,
        run=_normalize_task_run(forcer, function, bound),
        task=_tool_definition(forcer, function.path),
    )


def _emit_workflow_call(forcer, function, bound, satisfied=None):
    from swl.dag.tooldefs import _workflow_definition
    return _emit_step_call(
        forcer,
        function,
        bound,
        satisfied,
        run={},
        task=_workflow_definition(forcer, function),
        step_type='workflow',
    )


def _emit_mapped_step(forcer, fn, source, key=None):
    target = fn.function
    outputs = list(target.signature.outputs.keys())
    step_id = _step_id(forcer, target.name)
    bindings, map_info = _mapped_step_bindings(forcer, fn, source, key=key)
    from swl.dag.tooldefs import _forced_signature, _tool_definition, _workflow_definition
    signature = _forced_signature(forcer, fn)
    step = StepCall(
        id=step_id,
        path=target.path,
        source=source,
        outputs=outputs,
        bindings=bindings,
        task=_tool_definition(forcer, target.path) if target.kind == 'task' else _workflow_definition(forcer, target),
        deps=sorted(_step_dependencies(forcer, bindings or {'source': source})),
        type=target.kind,
        map=map_info,
        input_schema={name: normalize_swl_type(spec.type.value if getattr(spec, 'type', None) is not None else None) for name, spec in signature.inputs.items()} if signature is not None else None,
        output_schema={name: normalize_swl_type(spec.type.value if getattr(spec, 'type', None) is not None else None) for name, spec in signature.outputs.items()} if signature is not None else None,
    )
    forcer.steps.append(step)
    return step


def _map_ports(logical_source, input_names):
    if not isinstance(logical_source, TableSource):
        return list(input_names), []
    table_cols = set(logical_source.columns.keys())
    scatter = [name for name in input_names if name in table_cols]
    if scatter:
        return scatter, [name for name in input_names if name not in table_cols]
    return list(input_names), []


def _mapped_step_bindings(forcer, fn, source, key=None):
    logical_source = _logical_table_source(forcer, source)
    info = {'source': _binding_to_public_dict(logical_source)}
    if key is not None:
        info['group_by'] = key
    from swl.dag.tooldefs import _forced_signature
    signature = _forced_signature(forcer, fn)
    if signature is not None and signature.inputs:
        info['scatter'], info['broadcast'] = _map_ports(logical_source, sorted(signature.inputs.keys()))
    return {}, info


def _normalize_task_inputs(forcer, function, bound):
    names = list(function.signature.inputs.keys())
    if not names:
        return {}
    if not _can_project_record_fields(forcer, bound):
        return {names[0]: bound}
    return {name: _project_input(forcer, bound, name) for name in names}


def _normalize_task_run(forcer, function, bound):
    names = list(function.signature.run.keys())
    if not names or not _can_project_record_fields(forcer, bound):
        return {}
    result = {}
    for name in names:
        spec = function.signature.run[name]
        if _project_field(forcer, bound, name) is _SENTINEL:
            if getattr(spec, 'parsed_default', None) is None:
                raise ValueError(
                    f'Required run parameter "{name}" for task {function.name} '
                    f'({function.path}) has no value and no default'
                )
            continue
        result[name] = _project_input(forcer, bound, name)
    return result


def _step_id(forcer, name):
    count = forcer.step_name_counts.get(name, 0) + 1
    forcer.step_name_counts[name] = count
    if count == 1:
        return name
    return f'{name}_{count}'


def _step_call_key(forcer, function, inputs):
    return (
        function.path,
        tuple(sorted((name, _value_key(value)) for name, value in inputs.items())),
    )


def _validate_bindings(forcer, inputs, satisfied, function_name):
    if not satisfied:
        return
    bound_names = set(inputs.keys())
    missing = satisfied - bound_names
    if missing:
        raise ValueError(
            f'Step call to {function_name} has bindings {bound_names} '
            f'but checker predicted satisfied inputs: {satisfied}. '
            f'Missing: {missing}'
        )


def _step_dependencies(forcer, inputs):
    deps = set()
    for value in inputs.values():
        deps.update(_value_dependencies(forcer, value))
    return deps


def _walk_values(forcer, value):
    yield value
    if isinstance(value, Field):
        yield from _walk_values(forcer, value.source)
    elif isinstance(value, Merge):
        yield from _walk_values(forcer, value.left)
        yield from _walk_values(forcer, value.right)
    elif isinstance(value, Record):
        for item in value.fields.values():
            yield from _walk_values(forcer, item)


def _value_dependencies(forcer, value):
    deps = set()
    for item in _walk_values(forcer, value):
        if isinstance(item, Field) and isinstance(item.source, StepCall):
            deps.add(item.source.id)
    return deps


def _logical_table_source(forcer, source):
    if isinstance(source, Record):
        input_fields = {
            name: value
            for name, value in source.fields.items()
            if isinstance(value, Input)
        }
        if input_fields and len(input_fields) == len(source.fields):
            return TableSource('table', dict(input_fields))
    return source


def _binding_to_public_dict(value):
    return binding_to_dict(value)
