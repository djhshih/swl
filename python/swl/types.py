SWL_TYPE_ALIASES = {
    'File': 'file',
    'string': 'str',
}


def normalize_swl_type(value):
    if value is None:
        return None
    return SWL_TYPE_ALIASES.get(value, value)


def is_optional_type(value):
    return bool(value and isinstance(value, str) and value.endswith('?'))


def is_array_type(value):
    return bool(value and isinstance(value, str) and value.startswith('[') and value.endswith(']'))


def array_item_type(value):
    if is_array_type(value):
        return value[1:-1]
    return None


def to_array_type(value):
    if value is None or is_array_type(value):
        return value
    return f'[{value}]'


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
    optional = is_optional_type(value)
    if optional:
        value = value[:-1]
    if value == '[file]':
        base = {'type': 'array', 'items': 'File'}
        return ['null', base] if optional else base
    if value == '[str]':
        base = {'type': 'array', 'items': 'string'}
        return ['null', base] if optional else base
    if value == '[int]':
        base = {'type': 'array', 'items': 'int'}
        return ['null', base] if optional else base
    if value == '[float]':
        base = {'type': 'array', 'items': 'float'}
        return ['null', base] if optional else base
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
    }.get(swl_type, ('val', 'string'))
