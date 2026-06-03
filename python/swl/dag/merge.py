from swl.dag.node import Field, Input, Literal, Merge, Record, StepCall


def _canonicalize_merges(value):
    if isinstance(value, Merge):
        left = _canonicalize_merges(value.left)
        right = _canonicalize_merges(value.right)
        if isinstance(left, Record) and isinstance(right, Record):
            merged = dict(left.fields)
            merged.update(right.fields)
            return Record(merged)
        if isinstance(left, Record) and isinstance(right, Merge):
            right = _canonicalize_merges(right)
            if isinstance(right, Record):
                merged = dict(left.fields)
                merged.update(right.fields)
                return Record(merged)
            return Merge(left, right)
        if isinstance(right, Record) and isinstance(left, Merge):
            left = _canonicalize_merges(left)
            if isinstance(left, Record):
                merged = dict(left.fields)
                merged.update(right.fields)
                return Record(merged)
            return Merge(left, right)
        return Merge(left, right)
    if isinstance(value, Record):
        return Record({name: _canonicalize_merges(item) for name, item in value.fields.items()})
    if isinstance(value, Field):
        return Field(_canonicalize_merges(value.source), value.name)
    return value


def _flatten_value_terms(value):
    if isinstance(value, Merge):
        return _flatten_value_terms(value.left) + _flatten_value_terms(value.right)
    return [value]


def _normalize_output_value(value):
    if isinstance(value, Merge):
        terms = _flatten_value_terms(value)
        record_fields = {}
        others = []
        for term in terms:
            normalized = _normalize_output_value(term)
            if isinstance(normalized, Record):
                record_fields.update(normalized.fields)
            else:
                others.append(normalized)
        if not others:
            return Record(record_fields)
        result = Record(record_fields) if record_fields else others[0]
        if record_fields:
            start = 0
        else:
            start = 1
        for term in others[start:]:
            result = Merge(result, term)
        return result
    if isinstance(value, Record):
        return Record({name: _normalize_output_value(item) for name, item in value.fields.items()})
    if isinstance(value, Field):
        source = _normalize_output_value(value.source)
        if isinstance(source, Record) and len(source.fields) == 1:
            only = next(iter(source.fields.values()))
            if isinstance(only, Input):
                return Field(only, value.name)
        return Field(source, value.name)
    return value


def _value_key(value):
    if isinstance(value, Input):
        return ('input', value.name)
    if isinstance(value, Literal):
        return ('lit', value.value)
    if isinstance(value, Field):
        return ('field', _value_key(value.source), value.name)
    if isinstance(value, Merge):
        return _merge_key(value)
    if isinstance(value, Record):
        return ('record', _record_fields_key(value.fields))
    if isinstance(value, StepCall):
        return ('step', value.path, value.id)
    return ('other', type(value).__name__, repr(value))


def _merge_key(value):
    terms = _flatten_value_terms(value)
    record_fields = {}
    others = []
    for term in terms:
        if isinstance(term, Record):
            record_fields.update(term.fields)
        else:
            others.append(_value_key(term))
    parts = []
    if record_fields:
        parts.append(('record', _record_fields_key(record_fields)))
    parts.extend(others)
    if len(parts) == 1:
        return parts[0]
    return ('merge', tuple(parts))


def _record_fields_key(fields):
    return tuple(sorted((name, _value_key(value)) for name, value in fields.items()))
