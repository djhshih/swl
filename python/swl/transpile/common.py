import re as _re

from swl.dag.node import Field, Literal, OutputSpec


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


def classify_var(name, input_names, output_names, run_names):
    from swl.syntax.task.bash import _BUILTIN_VARS
    if name in _BUILTIN_VARS:
        return 'shell'
    if name in run_names:
        return 'run'
    if name in input_names:
        return 'input'
    if name in output_names:
        return 'output'
    return None



def interp_script(body, var_fn, expr_fn, joiner=''):
    pattern = r'(?<!\$)\$\{(.+?)\}|(?<!\$)\$(\w+)'
    result = []
    for line in body.split('\n'):
        if '$' not in line:
            result.append(line)
            continue

        def replace(m):
            if m.group(1) is not None:
                content = m.group(1).strip()
                if _re.fullmatch(r'\w+', content):
                    return var_fn(content) or m.group(0)
                return expr_fn(content) or m.group(0)
            return var_fn(m.group(2)) or m.group(0)

        result.append(_re.sub(pattern, replace, line))
    return '\n'.join(result)


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


def _flatten_output_names(outputs):
    """Replace dots with underscores in output names for flat targets (NF, SMK)."""
    return {
        name.replace('.', '_'): spec
        for name, spec in outputs.items()
    }


def format_resource_directives(run, formatters):
    directives = []
    for name, fmt in formatters.items():
        value = run_value(run, name)
        if value is not None:
            directives.append(fmt(value))
    return directives


def validate_dag_for_transpile(dag, backend_name):
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if isinstance(value, Literal):
            raise ValueError(f'{backend_name} does not support literal workflow outputs: {name}')


def flatten_dag_outputs(dag):
    dag.outputs = _flatten_output_names(dag.outputs)
    return dag


def field_chain_parts(value):
    current = value
    chain = []
    while isinstance(current, Field):
        chain.append(current.name)
        current = current.source
    return current, chain[-1] if chain else None, None


def field_path_after_first(value):
    _, _, tail = field_chain_parts(value)
    return tail or value.name



