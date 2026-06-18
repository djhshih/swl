import json
import re

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import (
    emit_name,
    flatten_dag_outputs,
    format_resource_directives,
    interp_script,
    run_value,
    step_name,
    validate_dag_for_transpile,
    validate_no_merge_bindings,
    word_interp,
)


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def transpile_dag_dict(data, workflow_id='main', _top_level=True, wrap_map=None):
    dag = DAG.from_dict(data)
    dag.validate()
    flatten_dag_outputs(dag)
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


def _output_path(step_id, leaf, scatter_ports):
    if scatter_ports:
        parts = [f'results/{step_id}']
        for port in sorted(scatter_ports):
            parts.append('{' + port + '}')
        parts.append(leaf)
        return '/'.join(parts)
    return f'results/{step_id}/{leaf}'


def _render_output_path(value, dag=None):
    if isinstance(value, Field) and isinstance(value.source, StepCall):
        prev = value.source
        out_name = value.name
        prev_scatter = set(prev.map_info.get('scatter', [])) if prev.is_mapped else set()
        spec = prev.task_def.get('outputs', {}).get(out_name, {})
        rendered = _output_template(spec)
        if rendered is not None:
            return _output_path(prev.id, rendered, prev_scatter)
        if prev.type == 'workflow':
            inner_default, inner_step_id = _inner_output_default(prev, out_name)
            if inner_default:
                rendered = _interp_to_smk(inner_default)
                if rendered is not None:
                    return _output_path(inner_step_id, rendered, set()) if inner_step_id else rendered
        if prev.is_mapped:
            scatter_ports = set(prev.map_info.get('scatter', []))
            return _output_path(prev.id, out_name, scatter_ports)
        return f'results/{prev.id}/{out_name}'
    return None



def _collect_wildcard_globs(dag):
    globs = []
    all_vars = set()
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if _is_mapped_output(value, dag):
            step = _step_for_output(value, dag)
            if step is not None:
                all_vars.update(_mapped_output_wildcards(step))
        else:
            rendered = _render_output_path(value, dag)
            if rendered:
                vars_in = re.findall(r'\{(\w+)\}', rendered)
                all_vars.update(vars_in)
    for var in sorted(all_vars):
        globs.append(f'{var.upper()} = config["{var}"]')
    return globs


def _inner_output_default(step, output_name):
    dag_data = step.task_def.get('dag', {})
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
                    return out_spec.get('path_template') or out_spec.get('default'), inner_step_id
    return None, None


def _mapped_output_wildcards(step):
    vars_set = set()
    dag_data = step.task_def.get('dag', {})
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
            scoped = _output_path(inner_step_id, rendered, set()) if inner_step_id else rendered
            vars_in = re.findall(r'\{(\w+)\}', scoped)
            if vars_in:
                kwargs = ', '.join(f'{v}={v.upper()}' for v in vars_in)
                return f'expand({repr(scoped)}, {kwargs})'
            return repr(scoped)
    scatter_ports = sorted(step.map_info.get('scatter', []))
    if not scatter_ports:
        return repr(f'results/{step.id}/{name}')
    first = scatter_ports[0]
    var = first.upper()
    dirs = [f'results/{step.id}', '{' + first + '}', name]
    pattern = '/'.join(dirs)
    return f'expand("{pattern}", {first}={var})'


def _validate_supported(dag):
    validate_dag_for_transpile(dag, 'Snakemake')
    validate_no_merge_bindings(dag, 'Snakemake')


def _rule_name(step_id):
    return step_name(step_id, 'rule')



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


def _output_template(spec):
    template = spec.get('path_template') if isinstance(spec, dict) else None
    if template is not None:
        return template
    default = spec.get('default') if isinstance(spec, dict) else None
    if default:
        return _interp_to_smk(default)
    return None


def _task_to_rule(step, dag, wrap_map=None):
    task = step.task_def
    body = task.get('body', '')
    rname = _rule_name(step.id)
    is_mapped = step.is_mapped

    scatter_ports = set()
    if is_mapped:
        group_by = step.map_info.get('group_by')
        if group_by:
            scatter_ports = {group_by}
        else:
            scatter_ports = set(step.map_info.get('scatter', []))

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
            path_expr = _binding_to_path(binding, in_name, is_mapped, scatter_ports, wrap_map, dag)
            lines.append(f'        {emit_name(in_name)}={path_expr},')
        lines.append('')

    task_outputs = task.get('outputs', {})
    params = _collect_params(step, dag)
    if task_outputs:
        lines.append('    output:')
        for out_name in task_outputs:
            rendered = _output_template(task_outputs[out_name])
            if rendered is not None:
                concrete = _resolve_concrete_path(rendered, params)
                if concrete is not None:
                    path_expr = repr(f'results/{step.id}/{concrete}')
                else:
                    path_expr = repr(_output_path(step.id, rendered, scatter_ports))
            else:
                path_expr = repr(_output_path(step.id, out_name, scatter_ports))
            lines.append(f'        {emit_name(out_name)}={path_expr},')
        lines.append('')

    params = _collect_params(step, dag)
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
        interp_body = _interpolate_shell(body, step, dag)
        for out_name, out_spec in task_outputs.items():
            pattern = _interp_preview(out_spec.get('default'))
            if pattern:
                m = re.match(r'\{(\w+)\}(.*)', pattern)
                if m:
                    var_name = m.group(1)
                    suffix = m.group(2)
                    param_ref = f'{{params.{var_name}}}{suffix}'
                    output_ref = f'{{output.{emit_name(out_name)}}}'
                    interp_body = interp_body.replace(param_ref, output_ref)
        lines.append('    shell:')
        lines.append('        """')
        for line in interp_body.split('\n'):
            lines.append(f'        {line}' if line.strip() else '')
        lines.append('        """')
        lines.append('')

    return '\n'.join(lines)



def _binding_to_path(binding, port_name, is_mapped, scatter_ports, wrap_map=None, dag=None):
    if binding is None:
        return f'config["{port_name}"]'

    if isinstance(binding, Input):
        return f'config["{binding.name}"]'

    if isinstance(binding, Field):
        if dag and binding.name in dag.inputs:
            return f'config["{binding.name}"]'
        if isinstance(binding.source, StepCall):
            prev_step = binding.source
            out_name = binding.name
            prev_scatter = set(prev_step.map_info.get('scatter', [])) if prev_step.is_mapped else set()
            if prev_step.type == 'workflow' and prev_step.is_mapped:
                inner_default, inner_step_id = _inner_output_default(prev_step, out_name)
                if inner_default:
                    rendered = _interp_to_smk(inner_default)
                    if rendered is not None:
                        scoped = _output_path(inner_step_id, rendered, set())
                        vars_in = re.findall(r'\{(\w+)\}', scoped)
                        if vars_in:
                            kwargs = ', '.join(f'{v}={v.upper()}' for v in vars_in)
                            return f'expand({repr(scoped)}, {kwargs})'
                        return repr(scoped)
            spec = prev_step.task_def.get('outputs', {}).get(out_name, {})
            rendered = _output_template(spec)
            if rendered is not None:
                return repr(_output_path(prev_step.id, rendered, prev_scatter))
            if prev_step.is_mapped:
                return repr(_output_path(prev_step.id, out_name, prev_scatter))
            return repr(f'results/{prev_step.id}/{out_name}')
        return f'config["{port_name}"]'

    if isinstance(binding, Literal):
        if isinstance(binding.value, str):
            return repr(binding.value)
        return json.dumps(binding.value)

    if isinstance(binding, Record):
        raise ValueError(
            f'Record binding with fields {list(binding.fields.keys())} must be flattened '
            'before Snakemake transpilation.'
        )

    if isinstance(binding, Merge):
        raise ValueError(
            f'Merge bindings must be flattened before Snakemake transpilation. '
            f'Found in binding {port_name}'
        )

    raise ValueError(f'Unsupported binding for Snakemake: {type(binding).__name__}')


def _interpolate_shell(body, step, dag=None):
    task = step.task_def
    task_inputs = task.get('inputs', {})
    input_names = set()
    for name, spec in task_inputs.items():
        if spec.get('type') not in ('str', 'int', 'float'):
            input_names.add(name)
    output_names = set(step.outputs)
    param_names = {p[0] for p in _collect_params(step, dag)}

    def _resolve_var(var):
        if var in input_names:
            return emit_name(var)
        if var in output_names:
            return emit_name(var)
        if var in param_names:
            return var
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


def _collect_params(step, dag=None):
    params = []
    task_inputs = step.task_def.get('inputs', {})
    for in_name, spec in task_inputs.items():
        typ = spec.get('type', 'str')
        if typ in ('str', 'int', 'float'):
            binding = step.bindings.get(in_name)
            if isinstance(binding, Literal):
                params.append((emit_name(in_name), json.dumps(binding.value)))
            elif isinstance(binding, Field) and isinstance(binding.source, StepCall):
                prev_step = binding.source
                out_name = binding.name
                prev_spec = prev_step.task_def.get('outputs', {}).get(out_name, {})
                default = prev_spec.get('default')
                if default:
                    rendered = _interp_to_smk(default)
                    params.append((emit_name(in_name), repr(rendered)) if rendered else (emit_name(in_name), repr(f'results/{prev_step.id}/{out_name}')))
                else:
                    params.append((emit_name(in_name), repr(f'results/{prev_step.id}/{out_name}')))
            else:
                params.append((emit_name(in_name), f'config["{in_name}"]'))
    run = step.task_def.get('run', {})
    for name in ('cpu', 'memory', 'time'):
        spec = run.get(name, {})
        if isinstance(spec, dict) and spec.get('type') in ('int', 'str', 'float', 'memory', 'time'):
            params.append((emit_name(name), json.dumps(spec['value'])))
    return params


def _emit_resources(step):
    run = step.task_def.get('run', {})
    directives = []

    directives.extend(format_resource_directives(run, {
        'cpu': lambda v: f'threads: {v}',
        'image': lambda v: f'container: "docker://{v}"',
    }))

    memory = run_value(run, 'memory')
    time_val = run_value(run, 'time')
    if memory is not None or time_val is not None:
        directives.append('resources:')
        if memory is not None:
            directives.append(f'    mem_mb={memory}')
        if time_val is not None:
            directives.append(f'    runtime_minutes={time_val}')

    return directives


def _dag_to_smk(dag, workflow_id, rule_names):
    lines = []

    mapped_outputs = {}
    simple_outputs = {}

    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        mapped = _is_mapped_output(value, dag)
        if mapped:
            mapped_outputs[name] = (value, output)
        else:
            simple_outputs[name] = (value, output)

    if simple_outputs or mapped_outputs:
        lines.append('rule all:')
        lines.append('    input:')

        for name, (value, output) in sorted(simple_outputs.items()):
            path_expr = _output_to_path(name, value, dag)
            lines.append(f'        {emit_name(name)}={path_expr},')

        for name, (value, output) in sorted(mapped_outputs.items()):
            step = _step_for_output(value, dag)
            if step is not None:
                lines.append(f'        {_mapped_output_expand(name, step)},')

        lines.append('')

    return '\n'.join(lines)


def _is_mapped_output(value, dag=None):
    if isinstance(value, Field) and isinstance(value.source, StepCall):
        return value.source.map is not None
    return False


def _step_for_output(value, dag=None):
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
            rendered = _render_output_path(value, dag)
            if rendered is not None:
                params = _collect_params(prev, dag)
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
    dag_data = step.task_def.get('dag', {})
    if not dag_data:
        return ''
    wf_id = f'{parent_id}_{step.id}'
    wrap_map = None
    if step.is_mapped:
        scatter_ports = step.map_info.get('scatter', [])
        if scatter_ports:
            wrap_map = {
                "step_id": step.id,
                "wildcard": scatter_ports[0],
                "scatter_ports": set(scatter_ports),
            }
    return transpile_dag_dict(dag_data, workflow_id=wf_id, _top_level=False, wrap_map=wrap_map)
