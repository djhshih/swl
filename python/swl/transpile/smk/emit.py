import json
import re

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import (
    emit_name,
    field_chain_parts,
    field_path_after_first,
    interp_script,
    run_value,
    source_kind,
    step_name,
    word_interp,
    workflow_name,
    _flatten_output_names,
)


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def transpile_dag_dict(data, workflow_id='main', _top_level=True, wrap_map=None):
    dag = DAG.from_dict(data)
    dag.validate()
    dag.outputs = _flatten_output_names(dag.outputs)
    _validate_supported(dag)

    rule_names = {}
    sub_workflows = []
    for step in dag.steps:
        if step.id not in rule_names:
            if step.type == 'workflow':
                sub = _subworkflow_to_smk(step, workflow_id)
                rule_names[step.id] = _rule_name(step.id)
                sub_workflows.append(sub)
            else:
                rule_names[step.id] = _rule_name(step.id)

    lines = []
    if _top_level:
        wildcard_globs = _collect_wildcard_globs(dag)
        if wildcard_globs:
            for g in wildcard_globs:
                lines.append(g)
            lines.append('')
        dag_lines = _dag_to_smk(dag, workflow_id, rule_names)
        if dag_lines.strip():
            lines.append(dag_lines)
            lines.append('')
    for step in dag.steps:
        if step.type == 'workflow':
            continue
        lines.append(_task_to_rule(step, dag, wrap_map))
        lines.append('')
    for sub in sub_workflows:
        lines.append(sub)
        lines.append('')
    return '\n'.join(lines)


def _resolve_concrete_path(pattern, params):
    vars_in = re.findall(r'\{(\w+)\}', pattern)
    if not vars_in:
        return pattern
    param_map = {name: expr for name, expr in params}
    replacements = {}
    for var in vars_in:
        expr = param_map.get(var)
        if expr is None:
            return None
        if not (expr.startswith('"') or expr.startswith("'")):
            return None
        replacements[var] = json.loads(expr)
    result = pattern
    for var, val in replacements.items():
        result = result.replace('{' + var + '}', str(val))
    return result


def _scope_path(step_id, rendered, is_mapped, scatter_ports):
    if is_mapped and scatter_ports:
        parts = [f'results/{step_id}']
        for port in sorted(scatter_ports):
            parts.append('{' + port + '}')
        parts.append(rendered)
        return repr('/'.join(parts))
    return repr(f'results/{step_id}/{rendered}')


def _scope_path_raw(step_id, rendered, is_mapped, scatter_ports):
    if is_mapped and scatter_ports:
        parts = [f'results/{step_id}']
        for port in sorted(scatter_ports):
            parts.append('{' + port + '}')
        parts.append(rendered)
        return '/'.join(parts)
    return f'results/{step_id}/{rendered}'


def _render_output_path(value):
    if isinstance(value, Field) and isinstance(value.source, StepCall):
        prev = value.source
        prev_is_mapped = prev.map is not None
        prev_scatter = set(prev.map.get('scatter', [])) if prev.map else set()
        spec = (prev.task or {}).get('outputs', {}).get(value.name, {})
        default_spec = spec.get('default')
        if default_spec:
            rendered = _interp_to_smk(default_spec)
            if rendered is not None:
                return _scope_path_raw(prev.id, rendered, prev_is_mapped, prev_scatter)
        if prev.type == 'workflow':
            inner_default, inner_step_id = _inner_output_default(prev, value.name)
            if inner_default:
                rendered = _interp_to_smk(inner_default)
                if rendered is not None:
                    return _scope_path_raw(inner_step_id, rendered, False, set()) if inner_step_id else rendered
        if prev.map is not None:
            scatter_ports = set(prev.map.get('scatter', []))
            return _default_output_path_raw(prev.id, value.name, scatter_ports)
        return f'results/{prev.id}/{value.name}'
    return None


def _default_output_path_raw(step_id, out_name, scatter_ports):
    if scatter_ports:
        parts = [f'results/{step_id}']
        for port in sorted(scatter_ports):
            parts.append('{' + port + '}')
        parts.append(out_name)
        return '/'.join(parts)
    return f'results/{step_id}/{out_name}'


def _collect_wildcard_globs(dag):
    globs = []
    all_vars = set()
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if _is_mapped_output(value):
            step = _step_for_output(value)
            if step is not None:
                all_vars.update(_mapped_output_wildcards(step))
        else:
            rendered = _render_output_path(value)
            if rendered:
                vars_in = re.findall(r'\{(\w+)\}', rendered)
                all_vars.update(vars_in)
    for var in sorted(all_vars):
        globs.append(f'{var.upper()} = config["{var}"]')
    return globs


def _inner_output_default(step, output_name):
    dag_data = (step.task or {}).get('dag', {})
    sub_outputs = dag_data.get('outputs', {})
    sub_value = sub_outputs.get(output_name, {})
    if isinstance(sub_value, dict):
        inner_value = sub_value.get('value', sub_value)
        inner_step_id = inner_value.get('step') if isinstance(inner_value, dict) else None
        inner_output_name = inner_value.get('output') if isinstance(inner_value, dict) else None
        if inner_step_id and inner_output_name:
            for s in dag_data.get('steps', []):
                if s.get('id') == inner_step_id:
                    out_spec = s.get('outputs', {}).get(inner_output_name, {})
                    return out_spec.get('default'), inner_step_id
    return None, None


def _mapped_output_wildcards(step):
    vars_set = set()
    dag_data = (step.task or {}).get('dag', {})
    for out_name in dag_data.get('outputs', {}):
        default_spec, _ = _inner_output_default(step, out_name)
        if default_spec:
            rendered = _interp_to_smk(default_spec)
            if rendered:
                vars_in = re.findall(r'\{(\w+)\}', rendered)
                vars_set.update(vars_in)
    return vars_set


def _mapped_output_expand(name, step):
    default_spec, inner_step_id = _inner_output_default(step, name)
    if default_spec:
        rendered = _interp_to_smk(default_spec)
        if rendered:
            scoped = _scope_path_raw(inner_step_id, rendered, False, set()) if inner_step_id else rendered
            vars_in = re.findall(r'\{(\w+)\}', scoped)
            if vars_in:
                kwargs = ', '.join(f'{v}={v.upper()}' for v in vars_in)
                return f'expand({repr(scoped)}, {kwargs})'
            return repr(scoped)
    scatter_ports = sorted(step.map.get('scatter', []))
    if not scatter_ports:
        return repr(f'results/{step.id}/{name}')
    first = scatter_ports[0]
    var = first.upper()
    dirs = [f'results/{step.id}', '{' + first + '}', name]
    pattern = '/'.join(dirs)
    return f'expand("{pattern}", {first}={var})'


def _validate_supported(dag):
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if isinstance(value, Literal):
            raise ValueError(f'Snakemake does not support literal workflow outputs: {name}')
    for step in dag.steps:
        for binding in step.bindings.values():
            if isinstance(binding, Merge):
                raise ValueError(
                    f'Merge bindings must be flattened before Snakemake transpilation. '
                    f'Found in step {step.id}'
                )


def _rule_name(step_id):
    return step_name(step_id, 'rule')


def _wf_name(workflow_id):
    return workflow_name(workflow_id, 'workflow')


def _san(name):
    return emit_name(name)


def _interp_preview(value):
    if isinstance(value, str):
        return value
    if isinstance(value, dict) and value.get('kind') == 'word':
        parts = value.get('parts', [])
        result = []
        for part in parts:
            kind = part.get('kind')
            if kind == 'literal':
                result.append(part.get('text', ''))
            elif kind == 'var':
                result.append('{' + part.get('name') + '}')
            elif kind == 'expr':
                result.append('${{' + part.get('text') + '}}')
        return ''.join(result)
    return None


def _interp_to_smk(value):
    if value is None:
        return None
    if isinstance(value, str):
        return value
    return word_interp(
        value,
        literal_fn=lambda text: text,
        var_fn=lambda name: '{' + name + '}',
        expr_fn=lambda text: '${{{' + text + '}}}',
    )


def _task_to_rule(step, dag, wrap_map=None):
    task = step.task or {}
    body = task.get('body', '')
    rname = _rule_name(step.id)
    is_mapped = step.map is not None

    scatter_ports = set()
    if is_mapped:
        group_by = step.map.get('group_by')
        if group_by:
            scatter_ports = {group_by}
        else:
            scatter_ports = set(step.map.get('scatter', []))

    lines = [f'rule {rname}:']

    task_inputs = task.get('inputs', {})
    file_inputs = {}
    string_inputs = {}
    for in_name, spec in task_inputs.items():
        typ = spec.get('type', 'file')
        if typ in ('str', 'int', 'float'):
            string_inputs[in_name] = spec
        else:
            file_inputs[in_name] = spec

    if file_inputs:
        lines.append('    input:')
        for in_name in file_inputs:
            binding = step.bindings.get(in_name)
            path_expr = _binding_to_path(binding, in_name, is_mapped, scatter_ports, wrap_map)
            lines.append(f'        {_san(in_name)}={path_expr},')
        lines.append('')

    task_outputs = task.get('outputs', {})
    params = _collect_params(step)
    if task_outputs:
        lines.append('    output:')
        for out_name in task_outputs:
            default_spec = task_outputs[out_name].get('default')
            if default_spec:
                rendered = _interp_to_smk(default_spec)
                if rendered is not None:
                    concrete = _resolve_concrete_path(rendered, params)
                    if concrete is not None:
                        path_expr = repr(f'results/{step.id}/{concrete}')
                    else:
                        path_expr = _scope_path(step.id, rendered, is_mapped, scatter_ports)
                else:
                    path_expr = _default_output_path(step.id, out_name, is_mapped, scatter_ports)
            else:
                path_expr = _default_output_path(step.id, out_name, is_mapped, scatter_ports)
            lines.append(f'        {_san(out_name)}={path_expr},')
        lines.append('')

    params = _collect_params(step)
    if params:
        lines.append('    params:')
        for name, expr in params:
            lines.append(f'        {name}={expr},')
        lines.append('')

    directives = _emit_resources(step)
    for d in directives:
        lines.append(f'    {d}')
    if directives:
        lines.append('')

    if body.strip():
        interp_body = _interpolate_shell(body, step)
        for out_name, out_spec in task_outputs.items():
            pattern = _interp_preview(out_spec.get('default'))
            if pattern:
                m = re.match(r'\{(\w+)\}(.*)', pattern)
                if m:
                    var_name = m.group(1)
                    suffix = m.group(2)
                    param_ref = f'{{params.{var_name}}}{suffix}'
                    output_ref = f'{{output.{_san(out_name)}}}'
                    interp_body = interp_body.replace(param_ref, output_ref)
        lines.append('    shell:')
        lines.append('        """')
        for line in interp_body.split('\n'):
            lines.append(f'        {line}' if line.strip() else '')
        lines.append('        """')
        lines.append('')

    return '\n'.join(lines)


def _default_output_path(step_id, out_name, is_mapped, scatter_ports):
    if is_mapped and scatter_ports:
        parts = [f'results/{step_id}']
        for port in sorted(scatter_ports):
            parts.append('{' + port + '}')
        parts.append(out_name)
        return repr('/'.join(parts))
    return repr(f'results/{step_id}/{out_name}')


def _binding_to_path(binding, port_name, is_mapped, scatter_ports, wrap_map=None):
    if binding is None:
        return f'config["{port_name}"]'

    if isinstance(binding, Input):
        return f'config["{binding.name}"]'

    if isinstance(binding, Field):
        if isinstance(binding.source, Input):
            root = binding.source
            return f'config["{root.name}"]'
        if isinstance(binding.source, StepCall):
            prev_step = binding.source
            prev_is_mapped = prev_step.map is not None
            prev_scatter = set(prev_step.map.get('scatter', [])) if prev_step.map else set()
            if prev_step.type == 'workflow' and prev_step.map is not None:
                inner_default, inner_step_id = _inner_output_default(prev_step, binding.name)
                if inner_default:
                    rendered = _interp_to_smk(inner_default)
                    if rendered is not None:
                        scoped = _scope_path_raw(inner_step_id, rendered, False, set())
                        vars_in = re.findall(r'\{(\w+)\}', scoped)
                        if vars_in:
                            kwargs = ', '.join(f'{v}={v.upper()}' for v in vars_in)
                            return f'expand({repr(scoped)}, {kwargs})'
                        return repr(scoped)
            spec = (prev_step.task or {}).get('outputs', {}).get(binding.name, {})
            default_spec = spec.get('default')
            if default_spec:
                rendered = _interp_to_smk(default_spec)
                if rendered is not None:
                    return _scope_path(prev_step.id, rendered, prev_is_mapped, prev_scatter)
            if prev_step.map is not None:
                return _default_output_path(prev_step.id, binding.name, True, prev_scatter)
            return repr(f'results/{prev_step.id}/{binding.name}')

    if isinstance(binding, Literal):
        if isinstance(binding.value, str):
            return repr(binding.value)
        return json.dumps(binding.value)

    if isinstance(binding, Record):
        raise ValueError(
            f'Record binding with fields {list(binding.fields.keys())} must be flattened '
            'before Snakemake transpilation.'
        )

    raise ValueError(f'Unsupported binding for Snakemake: {type(binding).__name__}')


def _interpolate_shell(body, step):
    task = step.task or {}
    task_inputs = task.get('inputs', {})
    input_names = set()
    for name, spec in task_inputs.items():
        if spec.get('type') not in ('str', 'int', 'float'):
            input_names.add(name)
    output_names = set(step.outputs)
    param_names = {p[0] for p in _collect_params(step)}

    def _resolve_var(var):
        if var in input_names:
            return '{input.' + _san(var) + '}'
        if var in output_names:
            return '{output.' + _san(var) + '}'
        if var in param_names:
            return '{params.' + var + '}'
        return None

    def var_fn(name):
        return _resolve_var(name)

    def expr_fn(text):
        resolved = text
        for var in re.findall(r'\w+', text):
            r = _resolve_var(var)
            if r:
                resolved = resolved.replace(var, r, 1)
        return '$(( ' + resolved + ' ))'

    return interp_script(body, var_fn, expr_fn, joiner='')


def _collect_params(step):
    params = []
    task_inputs = (step.task or {}).get('inputs', {})
    for in_name, spec in task_inputs.items():
        typ = spec.get('type', 'str')
        if typ in ('str', 'int', 'float'):
            binding = step.bindings.get(in_name)
            if isinstance(binding, Literal):
                params.append((_san(in_name), json.dumps(binding.value)))
            elif isinstance(binding, Field) and isinstance(binding.source, StepCall):
                prev_step = binding.source
                prev_spec = (prev_step.task or {}).get('outputs', {}).get(binding.name, {})
                default_spec = prev_spec.get('default')
                if default_spec:
                    rendered = _interp_to_smk(default_spec)
                    if rendered is not None:
                        params.append((_san(in_name), repr(rendered)))
                    else:
                        params.append((_san(in_name), repr(f'results/{prev_step.id}/{binding.name}')))
                else:
                    params.append((_san(in_name), repr(f'results/{prev_step.id}/{binding.name}')))
            else:
                params.append((_san(in_name), f'config["{in_name}"]'))
    run = (step.task or {}).get('run', {})
    for name in ('cpu', 'memory', 'time'):
        spec = run.get(name, {})
        if isinstance(spec, dict) and spec.get('type') in ('int', 'str', 'float', 'memory', 'time'):
            params.append((_san(name), json.dumps(spec['value'])))
    return params


def _emit_resources(step):
    run = (step.task or {}).get('run', {})
    directives = []

    cpu = run_value(run, 'cpu')
    if cpu is not None:
        directives.append(f'threads: {cpu}')

    memory = run_value(run, 'memory')
    time_val = run_value(run, 'time')
    has_resources = memory is not None or time_val is not None
    if has_resources:
        directives.append('resources:')
        if memory is not None:
            directives.append(f'    mem_mb={memory}')
        if time_val is not None:
            directives.append(f'    runtime_minutes={time_val}')

    image = run_value(run, 'image')
    if image is not None:
        directives.append(f'container: "docker://{image}"')

    return directives


def _dag_to_smk(dag, workflow_id, rule_names):
    lines = []

    mapped_outputs = {}
    simple_outputs = {}

    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        mapped = _is_mapped_output(value)
        if mapped:
            mapped_outputs[name] = (value, output)
        else:
            simple_outputs[name] = (value, output)

    if simple_outputs or mapped_outputs:
        lines.append('rule all:')
        lines.append('    input:')

        for name, (value, output) in sorted(simple_outputs.items()):
            path_expr = _output_to_path(name, value, dag)
            lines.append(f'        {_san(name)}={path_expr},')

        for name, (value, output) in sorted(mapped_outputs.items()):
            step = _step_for_output(value)
            if step is not None:
                lines.append(f'        {_mapped_output_expand(name, step)},')

        lines.append('')

    return '\n'.join(lines)


def _is_mapped_output(value):
    if isinstance(value, Field) and isinstance(value.source, StepCall):
        return value.source.map is not None
    return False


def _step_for_output(value):
    if isinstance(value, Field) and isinstance(value.source, StepCall):
        return value.source
    return None


def _rendered_to_fstring(rendered):
    parts = re.split(r'\{(\w+)\}', rendered)
    result = []
    for i, part in enumerate(parts):
        if i % 2 == 0:
            if part:
                result.append(repr(part))
        else:
            result.append(part.upper())
    if not result:
        return repr(rendered)
    if len(result) == 1:
        return result[0]
    return ' + '.join(result)


def _output_to_path(name, value, dag):
    if isinstance(value, Input):
        return f'config["{value.name}"]'
    if isinstance(value, Literal):
        return repr(value.value)
    if isinstance(value, Field):
        if isinstance(value.source, StepCall):
            prev = value.source
            rendered = _render_output_path(value)
            if rendered is not None:
                params = _collect_params(prev)
                concrete = _resolve_concrete_path(rendered, params) if params else None
                if concrete is not None:
                    return repr(concrete)
                vars_in = re.findall(r'\{(\w+)\}', rendered)
                if vars_in:
                    return _rendered_to_fstring(rendered)
                return repr(rendered)
            return repr(f'results/{prev.id}/{value.name}')
        if isinstance(value.source, Input):
            return f'config["{value.source.name}"]'
    if isinstance(value, Record):
        parts = []
        for fname, fval in sorted(value.fields.items()):
            parts.append(_output_to_path(f'{name}_{fname}', fval, dag))
        return ', '.join(parts)
    raise ValueError(f'Unsupported output value for Snakemake: {type(value).__name__}')


def _subworkflow_to_smk(step, parent_id):
    definition = step.task or {}
    dag_data = definition.get('dag', {})
    if not dag_data:
        return ''
    wf_id = f'{parent_id}_{step.id}'
    wrap_map = None
    if step.map is not None:
        scatter_ports = step.map.get('scatter', [])
        if scatter_ports:
            wrap_map = {
                "step_id": step.id,
                "wildcard": scatter_ports[0],
                "scatter_ports": set(scatter_ports),
            }
    return transpile_dag_dict(dag_data, workflow_id=wf_id, _top_level=False, wrap_map=wrap_map)
