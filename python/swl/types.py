SWL_TYPE_ALIASES = {
    'File': 'file',
    'string': 'str',
}


def normalize_swl_type(value):
    if value is None:
        return None
    if _is_invalid_array_of_optional(value):
        raise ValueError(
            f'Invalid type: {value!r}. '
            f'Use {value[:-1]}? for optional array, or {value[:-2]} for required array.'
        )
    return SWL_TYPE_ALIASES.get(value, value)


def _is_invalid_array_of_optional(value):
    return isinstance(value, str) and value.startswith('[') and value.endswith(']') and value[:-1].endswith('?')


def is_optional_type(value):
    return bool(value and isinstance(value, str) and value.endswith('?'))


def is_array_type(value):
    return bool(value and isinstance(value, str) and value.startswith('[') and value.endswith(']'))


def array_item_type(value):
    if is_array_type(value):
        return value[1:-1]
    return None


def to_array_type(value):
    if value is None:
        return None
    if _is_array_form(value):
        return value
    base = value[:-1] if value.endswith('?') else value
    return f'[{base}]'


def _is_array_form(value):
    if not isinstance(value, str) or not value.startswith('['):
        return False
    if value.endswith(']') and not value[:-1].endswith('?'):
        return True
    if value.endswith(']?') and not value[:-2].endswith('?'):
        return True
    return False


def base_scalar_type(value):
    if value is None:
        return None
    result = value
    if result.endswith('?'):
        result = result[:-1]
    if result.startswith('[') and result.endswith(']'):
        result = result[1:-1]
    return result


def to_cwl_type(value):
    if not value:
        return 'string'
    # Handle optional array: [type]?
    if value.endswith(']?'):
        inner = value[1:-2]
        item_cwl = to_cwl_type(inner)
        return ['null', {'type': 'array', 'items': item_cwl}]
    # Handle array: [type]
    if value.startswith('[') and value.endswith(']'):
        inner = value[1:-1]
        if inner.endswith('?'):
            raise ValueError(
                f'Invalid type: {value!r}. '
                f'[type?] is not allowed — use [type] for required array or [type]? for optional array.'
            )
        item_cwl = to_cwl_type(inner)
        return {'type': 'array', 'items': item_cwl}
    # Handle optional scalar
    optional = value.endswith('?')
    if optional:
        value = value[:-1]
    base = {
        'file': 'File',
        'str': 'string',
        'string': 'string',
        'int': 'int',
        'float': 'float',
        'bool': 'boolean',
        'boolean': 'boolean',
    }.get(value, 'string')
    return ['null', base] if optional else base


def to_wdl_type(swl_type, optional=False):
    if isinstance(swl_type, str) and swl_type.endswith('?'):
        optional = True
        swl_type = swl_type[:-1]
    base = {
        'file': 'File',
        'str': 'String',
        'int': 'Int',
        'float': 'Float',
        '[file]': 'Array[File]',
        '[str]': 'Array[String]',
        '[int]': 'Array[Int]',
        '[float]': 'Array[Float]',
    }.get(swl_type, 'String')
    if optional:
        return base + '?'
    return base


def to_nf_qualifier(swl_type):
    return {
        'file': ('path', None),
        'str': ('val', 'string'),
        'int': ('val', 'integer'),
        'float': ('val', 'float'),
        '[file]': ('path', None),
    }.get(swl_type, ('val', 'string'))
