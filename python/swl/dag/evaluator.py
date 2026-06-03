import swl.ir.node as ir
from swl.dag.context import ForceEnv, _SENTINEL
from swl.dag.merge import _canonicalize_merges, _value_key
from swl.dag.node import Field, ForcedFunction, Input, Literal, Merge, Record, StepCall


def _forced_function_key(value):
    function = value.function
    if isinstance(function, ir.Function):
        function_key = ('function', function.kind, function.path)
    elif isinstance(function, ir.Lambda):
        function_key = ('lambda', id(function))
    else:
        return None
    return (function_key, _value_key(value.bound))


def _assert_canonical(forcer, node):
    children = _canonical_children(node)
    if children is None:
        return
    for child in children:
        _assert_canonical(forcer, child)



def _canonical_children(node):
    if isinstance(node, ir.Lambda):
        return [node.body]
    if isinstance(node, ir.Block):
        return [*node.bindings, node.result]
    if isinstance(node, ir.Variable):
        return [node.value]
    if isinstance(node, (ir.Apply, ir.Map, ir.Update)):
        return [node.function, node.arg] if hasattr(node, 'function') else [node.left, node.right]
    if isinstance(node, ir.Record):
        return list(node.fields.values())
    if isinstance(node, ir.Field):
        return [node.record]
    return None


def force_value(forcer, node, env):
    if isinstance(node, ir.Literal):
        return Literal(node.value)

    if isinstance(node, ir.Input):
        if node.name not in forcer.inputs:
            optional = bool(node.type and node.type.endswith('?'))
            forcer.inputs[node.name] = Input(node.name, node.type, node.desc, optional)
        return forcer.inputs[node.name]

    if isinstance(node, ir.Name):
        value = env.lookup(node.name)
        if value is not None:
            return value
        raise ValueError(f'Undefined variable during forcing: {node.name}')

    if isinstance(node, ir.Ref):
        return _force_ref(forcer, node, env)

    if isinstance(node, ir.Record):
        return Record({name: force_value(forcer, value, env) for name, value in node.fields.items()})

    if isinstance(node, ir.Field):
        source = force_value(forcer, node.record, env)
        projected = _project_field(forcer, source, node.name)
        if projected is not _SENTINEL:
            return projected
        if isinstance(source, StepCall):
            if node.name not in source.outputs:
                raise ValueError(f'Missing field on tab: {node.name}')
            return Field(source, node.name)
        return Field(source, node.name)

    if isinstance(node, ir.Update):
        left = force_value(forcer, node.left, env)
        right = force_value(forcer, node.right, env)
        if isinstance(left, StepCall):
            raise ValueError(
                f'Record update (//) on a task/workflow call result is not supported: '
                f'left operand is a step call ({left.id})'
            )
        if isinstance(right, StepCall):
            raise ValueError(
                f'Record update (//) on a task/workflow call result is not supported: '
                f'right operand is a step call ({right.id})'
            )
        if isinstance(left, Input) and isinstance(right, Record):
            raise ValueError(
                'Table-record update (tab // rec) is not supported during forcing: '
                f'a table input ({left.name}) cannot be merged with a record literal. '
                'Use explicit bindings instead of // on table inputs.'
            )
        if isinstance(left, Record) and isinstance(right, Input):
            raise ValueError(
                'Record-table update (rec // tab) is not supported during forcing: '
                f'a table input ({right.name}) cannot be merged with a record literal. '
                'Use explicit bindings instead of // on table inputs.'
            )
        return _merge_values(forcer, left, right)

    if isinstance(node, ir.Function):
        return ForcedFunction(node, None, node.signature)

    if isinstance(node, ir.Lambda):
        return ForcedFunction(node, None, node.signature)

    if isinstance(node, ir.Apply):
        return _force_apply(forcer, node, env)

    if isinstance(node, ir.Map):
        fn = force_value(forcer, node.function, env)
        arg = force_value(forcer, node.arg, env)
        return _force_map(forcer, fn, arg, key=getattr(node, 'key', None))

    if isinstance(node, ir.Block):
        local = ForceEnv(env)
        _register_block(forcer, node.bindings)
        for bind in node.bindings:
            local.bind(bind.name, ir.Ref(bind.id, bind.name))
        return force_value(forcer, node.result, local)

    if isinstance(node, ir.Variable):
        forcer.variables[node.id] = node
        value = force_value(forcer, node.value, env)
        forcer.forced_variables[node.id] = value
        return value

    if isinstance(node, ir.Unknown):
        raise ValueError(f'Unknown IR node reached forcing (the lowerer should have rejected this): {node!r}')

    raise ValueError(f'Unsupported IR node during forcing: {type(node).__name__}')


def _force_ref(forcer, ref, env):
    if ref.id in forcer.forced_variables:
        return forcer.forced_variables[ref.id]
    variable = forcer.variables[ref.id]
    value = force_value(forcer, variable.value, env)
    forcer.forced_variables[ref.id] = value
    return value


def _force_apply(forcer, node, env):
    from swl.dag.emit import _emit_mapped_step
    fn = force_value(forcer, node.function, env)
    arg = force_value(forcer, node.arg, env)
    key = _apply_key(forcer, node.function, fn, arg)
    if key is not None and key in forcer.apply_cache:
        return forcer.apply_cache[key]
    satisfied = getattr(node, 'satisfied', set())
    result = _apply(forcer, fn, arg, satisfied)
    result = _canonicalize_merges(result)
    if key is not None:
        forcer.apply_cache[key] = result
    return result


def _apply(forcer, fn, arg, satisfied=None):
    if satisfied is None:
        satisfied = set()
    builtin = _apply_builtin_map(forcer, fn, arg, satisfied)
    if builtin is not _SENTINEL:
        return builtin
    if isinstance(arg, StepCall):
        return _apply_mapped(forcer, fn, arg)
    if not isinstance(fn, ForcedFunction):
        raise ValueError(f'Cannot apply non-function value during forcing: {fn!r}')

    bound = _merge_bound(forcer, fn.bound, arg)
    accumulated = fn.satisfied | satisfied

    if isinstance(fn.function, ir.Function):
        emitted = _apply_emittable_function(forcer, fn.function, bound, accumulated)
        if emitted is not _SENTINEL:
            return emitted

    if isinstance(fn.function, ir.Lambda):
        if not isinstance(bound, Record):
            return ForcedFunction(fn.function, bound, fn.signature, accumulated)
        return _force_lambda(forcer, fn.function, bound)

    return ForcedFunction(fn.function, bound, fn.signature, accumulated)



def _apply_builtin_map(forcer, fn, arg, satisfied):
    if not (isinstance(fn, ForcedFunction) and isinstance(fn.function, ir.Function)):
        return _SENTINEL
    if fn.function.kind != 'builtin' or fn.function.name not in ('map', 'map_by'):
        return _SENTINEL
    if fn.function.name == 'map':
        if fn.bound is None:
            return ForcedFunction(fn.function, Record({'f': arg}), fn.signature, satisfied)
        if isinstance(fn.bound, Record) and 'f' in fn.bound.fields:
            return _force_map(forcer, fn.bound.fields['f'], arg)
        return ForcedFunction(fn.function, _merge_bound(forcer, fn.bound, arg), fn.signature, satisfied)
    if fn.bound is None:
        return ForcedFunction(fn.function, Record({'f': arg}), fn.signature, satisfied)
    if isinstance(fn.bound, Record) and 'f' in fn.bound.fields and 'key' not in fn.bound.fields:
        return ForcedFunction(fn.function, Record({'f': fn.bound.fields['f'], 'key': arg}), fn.signature, satisfied)
    if isinstance(fn.bound, Record) and 'f' in fn.bound.fields and 'key' in fn.bound.fields:
        key_value = fn.bound.fields['key']
        key = key_value.value if isinstance(key_value, Literal) else None
        return _force_map(forcer, fn.bound.fields['f'], arg, key=key)
    return ForcedFunction(fn.function, _merge_bound(forcer, fn.bound, arg), fn.signature, satisfied)



def _apply_emittable_function(forcer, function, bound, accumulated):
    if function.kind == 'task':
        if not _saturated_task(forcer, function, bound):
            return _SENTINEL
        from swl.dag.emit import _emit_task_call
        return _emit_task_call(forcer, function, bound, accumulated)
    if function.kind == 'workflow':
        if not _saturated_signature(forcer, function, bound):
            return _SENTINEL
        from swl.dag.emit import _emit_workflow_call
        return _emit_workflow_call(forcer, function, bound, accumulated)
    return _SENTINEL


def _force_lambda(forcer, function, bound):
    local = ForceEnv()
    local.bind(function.param, bound)
    return force_value(forcer, function.body, local)


def _force_map(forcer, fn, arg, key=None):
    if not isinstance(fn, ForcedFunction):
        raise ValueError(f'Cannot map non-function value during forcing: {fn!r}')
    source = arg
    if isinstance(arg, Record) and len(arg.fields) == 1:
        source = next(iter(arg.fields.values()))
    return _apply_mapped(forcer, fn, source, key=key)


def _apply_mapped(forcer, fn, source, key=None):
    if not isinstance(fn, ForcedFunction):
        raise ValueError(f'Cannot map non-function value during forcing: {fn!r}')
    if _is_batch_function(forcer, fn):
        raise ValueError('map on batch workflow is not supported during forcing')
    if isinstance(fn.function, ir.Function) and fn.function.kind in ('task', 'workflow'):
        from swl.dag.emit import _emit_mapped_step
        return _emit_mapped_step(forcer, fn, source, key=key)
    raise ValueError(f'map requires normalized executable callable during forcing: {fn.function!r}')


def _is_batch_function(forcer, fn):
    function = fn.function if isinstance(fn, ForcedFunction) else fn
    return bool(getattr(function, 'is_batch', False))


def _register_block(forcer, bindings):
    for bind in bindings:
        forcer.variables[bind.id] = bind


def _merge_bound(forcer, old, new):
    if old is None:
        return new
    return _merge_values(forcer, old, new)


def _available_inputs(forcer, value):
    if isinstance(value, Input):
        return {value.name}
    if isinstance(value, Record):
        return set(value.fields.keys())
    if isinstance(value, Merge):
        return _available_inputs(forcer, value.left).union(_available_inputs(forcer, value.right))
    if isinstance(value, StepCall):
        return set(value.outputs)
    return set()


def _is_saturated(forcer, function, bound):
    if _is_opaque_record_carrier(forcer, bound):
        return True
    available = _available_inputs(forcer, bound)
    return all(name in available for name in function.signature.inputs.keys())


def _saturated_signature(forcer, function, bound):
    return _is_saturated(forcer, function, bound)


def _saturated_task(forcer, function, bound):
    return _is_saturated(forcer, function, bound)


def _forced_apply_key(forcer, value):
    if not isinstance(value, ForcedFunction):
        return None
    return _forced_function_key(value)


def _apply_key(forcer, function_node, fn, arg):
    arg_key = _value_key(arg)
    forced_key = _forced_apply_key(forcer, fn)
    if forced_key is not None:
        return ('apply-forced', forced_key, arg_key)
    if isinstance(function_node, ir.Ref):
        return ('apply-ref', function_node.id, arg_key)
    return None


def _looks_record_like(forcer, value):
    return isinstance(value, (Record, Merge))


def _can_project_record_fields(forcer, value):
    return isinstance(value, (Input, Field, Record, Merge))


def _is_opaque_record_carrier(forcer, value):
    if isinstance(value, (Input, Field)):
        return True
    if isinstance(value, Record) and len(value.fields) == 1:
        only = next(iter(value.fields.values()))
        return isinstance(only, (Input, Field))
    return False


def _project_field(forcer, source, name):
    if isinstance(source, (Input, Field)):
        return Field(source, name)
    if isinstance(source, Record):
        if name in source.fields:
            return source.fields[name]
        if len(source.fields) == 1:
            only = next(iter(source.fields.values()))
            if isinstance(only, (Input, Field)):
                return Field(only, name)
        return _SENTINEL
    if isinstance(source, Merge):
        projected = _project_field(forcer, source.right, name)
        if projected is not _SENTINEL:
            return projected
        return _project_field(forcer, source.left, name)
    return _SENTINEL


def _merge_values(forcer, left, right):
    if isinstance(left, Record) and isinstance(right, Record):
        fields = dict(left.fields)
        fields.update(right.fields)
        return Record(fields)
    return Merge(left, right)


def _project_input(forcer, value, name):
    if _can_project_record_fields(forcer, value):
        projected = _project_field(forcer, value, name)
        if projected is not _SENTINEL:
            return projected
        return Field(value, name)
    return value
