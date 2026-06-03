from swl.dag.merge import _canonicalize_merges, _flatten_value_terms, _normalize_output_value, _value_key
from swl.dag.node import DAG, Field, ForcedFunction, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.dag.tooldefs import _forced_signature, _input
from swl.dag.emit import _walk_values
from swl.types import is_optional_type, to_array_type


def _finalize_dag(forcer, value):
    _refine_input_metadata(forcer)
    _flatten_step_bindings(forcer)
    outputs = _final_outputs(forcer, value)
    outputs = _flatten_outputs(forcer, outputs)
    output_specs = _build_output_specs(forcer, outputs)
    for name, output in output_specs.items():
        _assert_wireable_output(forcer, name, output.value)
    _prune_unused_inputs(forcer, value, output_specs)
    dag = DAG(dict(forcer.inputs), list(forcer.steps), output_specs)
    dag.validate()
    return dag


def _refine_input_metadata(forcer):
    refined = {}
    for name, current in forcer.inputs.items():
        candidates = []
        for step in forcer.steps:
            spec = step.task.get('inputs', {}).get(name)
            if spec is not None:
                candidates.append(spec)
        best = current
        if candidates:
            typ = _merge_input_type(forcer, name, current.type, candidates)
            desc = _merge_input_desc(forcer, current.desc, candidates)
            optional = bool(typ and typ.endswith('?'))
            best = Input(name, typ, desc, optional)
        refined[name] = best
    forcer.inputs = refined


def _merge_input_type(forcer, name, current, candidates):
    types = []
    if current is not None:
        types.append(current)
    types.extend(spec.get('type') for spec in candidates if spec.get('type') is not None)
    unique = list(dict.fromkeys(types))
    if not unique:
        return None
    if len(unique) > 1:
        array_unique = list(dict.fromkeys(to_array_type(value) for value in unique))
        if len(array_unique) == 1:
            return array_unique[0]
        raise ValueError(f'Conflicting input types during forcing: {name}: {unique}')
    return unique[0]


def _merge_input_desc(forcer, current, candidates):
    descs = []
    if current:
        descs.append(current)
    descs.extend(spec.get('desc') for spec in candidates if spec.get('desc'))
    unique = list(dict.fromkeys(descs))
    if not unique:
        return None
    return ' / '.join(unique)


def _flatten_step_bindings(forcer):
    for step in forcer.steps:
        flat = {}
        for name, value in step.bindings.items():
            flattened = _flatten_merge_value(forcer, value)
            if isinstance(flattened, Record):
                flat.update(flattened.fields)
            else:
                flat[name] = flattened
        step.bindings = flat


def _flatten_outputs(forcer, outputs):
    flat = {}
    for name, value in outputs.items():
        flattened = _flatten_merge_value(forcer, value)
        if isinstance(flattened, Record):
            flat.update(flattened.fields)
        else:
            flat[name] = flattened
    return flat


def _build_output_specs(forcer, outputs):
    specs = {}
    for name, value in outputs.items():
        output_type = _infer_output_type(forcer, value)
        desc = None
        normalized = _normalize_output_value(value)
        if isinstance(normalized, Field) and isinstance(normalized.source, StepCall):
            spec = (normalized.source.task or {}).get('outputs', {}).get(normalized.name, {})
            desc = spec.get('desc')
        specs[name] = OutputSpec(
            type=output_type,
            desc=desc,
            optional=is_optional_type(output_type) if output_type else False,
            value=value,
        )
    return specs


def _infer_output_type(forcer, value):
    normalized = _normalize_output_value(value)
    if isinstance(normalized, Input):
        return normalized.type or 'str'
    if isinstance(normalized, Literal):
        if isinstance(normalized.value, bool):
            return 'bool'
        if isinstance(normalized.value, int):
            return 'int'
        if isinstance(normalized.value, float):
            return 'float'
        if isinstance(normalized.value, str):
            return 'str'
        return 'str'
    if isinstance(normalized, Record):
        return None
    if isinstance(normalized, Field):
        if isinstance(normalized.source, Input):
            source_type = normalized.source.type
            if source_type and source_type.startswith('[') and source_type.endswith(']'):
                return source_type[1:-1]
            return source_type or 'str'
        if isinstance(normalized.source, StepCall):
            spec = (normalized.source.task or {}).get('outputs', {}).get(normalized.name, {})
            typ = spec.get('type')
            if getattr(normalized.source, 'map', None) is not None and typ is not None:
                return to_array_type(typ)
            return typ or 'str'
        if isinstance(normalized.source, Field):
            current = normalized
            while isinstance(current, Field) and isinstance(current.source, Field):
                current = current.source
            if isinstance(current, Field) and isinstance(current.source, StepCall):
                spec = (current.source.task or {}).get('outputs', {}).get(current.name, {})
                typ = spec.get('type')
                if getattr(current.source, 'map', None) is not None and typ is not None:
                    return to_array_type(typ)
                return typ or 'str'
            return 'str'
    return 'str'


def _flatten_merge_value(forcer, value):
    if isinstance(value, Merge):
        terms = _flatten_value_terms(value)
        record_fields = {}
        non_record = []
        for term in terms:
            flattened = _flatten_merge_value(forcer, term)
            if isinstance(flattened, Record):
                record_fields.update(flattened.fields)
            else:
                non_record.append(flattened)
        if not non_record:
            return Record(record_fields)
        if len(non_record) == 1 and not record_fields:
            return non_record[0]
        if len(non_record) == 1 and record_fields:
            return Record({**record_fields, **_unwrap_non_record(forcer, non_record[0])})
        if not record_fields:
            raise ValueError(
                f'Cannot flatten merge: incompatible non-record values '
                f'({len(non_record)} terms)'
            )
        return Record(record_fields)
    if isinstance(value, Record):
        return Record({
            name: _flatten_merge_value(forcer, item)
            for name, item in value.fields.items()
        })
    return value


def _unwrap_non_record(forcer, value):
    if isinstance(value, Input):
        return {}
    if isinstance(value, Field):
        return {value.name: value}
    return {}


def _prune_unused_inputs(forcer, root_value, outputs):
    used = set()
    if isinstance(root_value, ForcedFunction):
        signature = _forced_signature(forcer, root_value)
        if signature is not None:
            used.update(signature.inputs.keys())
    for step in forcer.steps:
        for value in step.bindings.values():
            for item in _walk_values(forcer, value):
                if isinstance(item, Input):
                    used.add(item.name)
        mapped_source = getattr(step, 'source', None)
        if mapped_source is not None:
            mapped_ports = getattr(step, 'map', {}).get('ports') if getattr(step, 'map', None) is not None else []
            for item in _walk_values(forcer, mapped_source):
                if isinstance(item, Input) and not mapped_ports:
                    used.add(item.name)
    for output in outputs.values():
        for item in _walk_values(forcer, output.value if isinstance(output, OutputSpec) else output):
            if isinstance(item, Input):
                used.add(item.name)
    for step in forcer.steps:
        mapped = getattr(step, 'map', None)
        if mapped is not None:
            source = mapped.get('source', {})
            ports = mapped.get('ports') or []
            if source.get('source') == 'input' and source.get('name') not in forcer.inputs and not ports:
                _input(forcer, source['name'])
            if source.get('source') == 'input' and not ports:
                used.add(source['name'])
    forcer.inputs = {name: value for name, value in forcer.inputs.items() if name in used}


def _assert_wireable_output(forcer, name, value):
    normalized = _normalize_output_value(value)
    if isinstance(normalized, (Input, Literal)):
        return
    if isinstance(normalized, Field) and isinstance(normalized.source, (Input, StepCall)):
        return
    raise ValueError(f'Workflow output did not normalize to a wireable value: {name}: {normalized!r}')


def _final_outputs(forcer, value):
    value = _canonicalize_merges(value)
    fields = _collect_output_fields(forcer, value)
    if fields is not None:
        return fields
    return {'result': value}


def _collect_output_fields(forcer, value):
    normalized = _normalize_output_value(value)
    if isinstance(normalized, Record):
        return dict(normalized.fields)
    return None
