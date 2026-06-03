from swl.dag.node import Field


def normalize_identifier(name, default, *, upper=False):
    value = name.replace('-', '_').lstrip('_')
    if not value:
        value = default
    if value[0].isdigit():
        value = '_' + value
    return value.upper() if upper else value.lower()


def workflow_name(name, default='main', *, upper=False):
    return normalize_identifier(name, default, upper=upper)


def step_name(name, default='task', *, upper=False):
    return normalize_identifier(name, default, upper=upper)


def emit_name(name):
    return name.replace('-', '_')


def field_chain_parts(value):
    chain = []
    current = value
    while isinstance(current, Field):
        chain.append(current)
        current = current.source
    if not chain:
        return value, None, None
    chain.reverse()
    return chain[0].source, chain[0].name, '.'.join(part.name for part in chain[1:])


def field_path_after_first(value):
    if not isinstance(value, Field):
        return None
    _, _, tail = field_chain_parts(value)
    return tail or None


def word_interp(value, literal_fn, var_fn, expr_fn, joiner=''):
    if value is None:
        return None
    if value.get('kind') != 'word':
        return None
    result = []
    for part in value.get('parts', []):
        kind = part.get('kind')
        if kind == 'literal':
            result.append(literal_fn(part.get('text', '')))
        elif kind == 'var':
            result.append(var_fn(part.get('name')))
        elif kind == 'expr':
            result.append(expr_fn(part.get('text')))
        else:
            raise ValueError(f'Unsupported interpolation part: {part!r}')
    return joiner.join(result)


def run_value(run, name):
    spec = run.get(name, {})
    if not isinstance(spec, dict):
        return None
    return spec.get('value')


def source_kind(source):
    return source.get('source') if isinstance(source, dict) else None


def table_columns(source):
    if source_kind(source) == 'table':
        return source.get('columns', {})
    return {}


def source_input_name(source):
    if source_kind(source) == 'input':
        return source.get('name')
    return None


def column_input_name(source, col_name):
    if source_kind(source) == 'input':
        return source['name']
    if source_kind(source) == 'table':
        col = source.get('columns', {}).get(col_name, {})
        if isinstance(col, dict) and col.get('source') == 'input':
            return col['name']
    return col_name
