import json
import re

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import (
    column_input_name,
    emit_name,
    field_chain_parts,
    field_path_after_first,
    run_value,
    source_input_name,
    source_kind,
    step_name,
    table_columns,
    word_interp,
    workflow_name,
)


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def transpile_dag_dict(data, workflow_id='main', _top_level=True):
    dag = DAG.from_dict(data)
    dag.validate()
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
        wildcard_lists = _collect_wildcard_lists(dag)
        if wildcard_lists:
            lines.append('# TODO: Define wildcard lists for mapped outputs')
            for port in sorted(wildcard_lists):
                lines.append(f'{port.upper()}_LIST = [...]')
            lines.append('')
    for step in dag.steps:
        if step.type == 'workflow':
            continue
        lines.append(_task_to_rule(step, dag))
        lines.append('')
    for sub in sub_workflows:
        lines.append(sub)
        lines.append('')
    if _top_level:
        lines.append(_dag_to_smk(dag, workflow_id, rule_names))
    return '\n'.join(lines)


def _collect_wildcard_lists(dag):
    ports = set()
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if _is_mapped_output(value):
            step = _step_for_output(value)
            if step is not None:
                for port in step.map.get('scatter', []):
                    ports.add(port)
    return ports


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


def _interp_to_smk(value):
    if value is None:
        return None
    return word_interp(
        value,
        literal_fn=lambda text: text,
        var_fn=lambda name: '{' + name + '}',
        expr_fn=lambda text: '${{{' + text + '}}}',
    )


def _task_to_rule(step, dag):
    task = step.task or {}
    body = task.get('body', '')
    rname = _rule_name(step.id)
    is_mapped = step.map is not None
    is_map_by = is_mapped and step.map.get('group_by') is not None

    if is_map_by:
        return _mapped_by_step_to_smk(step, dag)

    scatter_ports = set()
    if is_mapped:
        scatter_ports = set(step.map.get('scatter', []))

    lines = [f'rule {rname}:', '']

    task_inputs = task.get('inputs', {})
    if task_inputs:
        lines.append('    input:')
        for in_name in task_inputs:
            binding = step.bindings.get(in_name)
            if is_mapped and in_name in scatter_ports:
                path_expr = repr('inputs/{' + in_name + '}')
            else:
                path_expr = _binding_to_path(binding, in_name, is_mapped, scatter_ports)
            lines.append(f'        {_san(in_name)}={path_expr},')
        lines.append('')

    task_outputs = task.get('outputs', {})
    if task_outputs:
        lines.append('    output:')
        for out_name in task_outputs:
            default_spec = task_outputs[out_name].get('default')
            if default_spec:
                rendered = _interp_to_smk(default_spec)
                if rendered is not None:
                    path_expr = repr(rendered)
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
        lines.append('    shell:')
        if '\n' in interp_body:
            lines.append('        """')
            for line in interp_body.split('\n'):
                lines.append(f'        {line}' if line.strip() else '')
            lines.append('        """')
        else:
            lines.append(f'        "{interp_body}"')
        lines.append('')

    lines.append('}')
    return '\n'.join(lines)


def _default_output_path(step_id, out_name, is_mapped, scatter_ports):
    if is_mapped and scatter_ports:
        parts = [f'results/{step_id}']
        for port in sorted(scatter_ports):
            parts.append('{' + port + '}')
        parts.append(out_name)
        return repr('/'.join(parts))
    return repr(f'results/{step_id}/{out_name}')


def _binding_to_path(binding, port_name, is_mapped, scatter_ports):
    if binding is None:
        if is_mapped and port_name in scatter_ports:
            return repr('inputs/{' + port_name + '}')
        return repr(f'inputs/{port_name}')

    if isinstance(binding, Input):
        if is_mapped and port_name in scatter_ports:
            return repr('inputs/{' + port_name + '}')
        return repr(f'inputs/{binding.name}')

    if isinstance(binding, Field):
        if isinstance(binding.source, Input):
            root = binding.source
            tail = field_path_after_first(binding)
            if tail:
                return repr(f'inputs/{root.name}/{tail}')
            return f'inputs.{root.name}'
        if isinstance(binding.source, StepCall):
            prev_step = binding.source
            spec = (prev_step.task or {}).get('outputs', {}).get(binding.name, {})
            default_spec = spec.get('default')
            if default_spec:
                rendered = _interp_to_smk(default_spec)
                if rendered is not None:
                    return repr(rendered)
            if prev_step.map is not None:
                prev_scatter = set(prev_step.map.get('scatter', []))
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
    input_names = set(task.get('inputs', {}).keys())
    output_names = set(step.outputs)
    param_names = {p[0] for p in _collect_params(step)}

    def replace_var(m):
        var = m.group(1) or m.group(2)
        if var in input_names:
            return '{input.' + _san(var) + '}'
        if var in output_names:
            return '{output.' + _san(var) + '}'
        if var in param_names:
            return '{params.' + var + '}'
        return m.group(0)

    pattern = r'(?<!\$)\$\{(\w+)\}|(?<!\$)\$(\w+)'
    result = []
    for line in body.split('\n'):
        result.append(re.sub(pattern, replace_var, line))
    return '\n'.join(result)


def _collect_params(step):
    params = []
    task_inputs = (step.task or {}).get('inputs', {})
    for in_name, spec in task_inputs.items():
        typ = spec.get('type', 'str')
        if typ in ('str', 'int', 'float'):
            binding = step.bindings.get(in_name)
            if isinstance(binding, Literal):
                params.append((_san(in_name), json.dumps(binding.value)))
            elif isinstance(binding, Input):
                params.append((_san(in_name), f'config["{binding.name}"]'))
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
                scatter_ports = sorted(step.map.get('scatter', []))
                dirs = [f'results/{step.id}']
                for port in scatter_ports:
                    dirs.append('{' + port + '}')
                dirs.append(name)
                pattern = '/'.join(dirs)
                expand_args = ', '.join(
                    f'{port}={port.upper() + "_LIST"}' for port in scatter_ports
                )
                lines.append(f'        expand("{pattern}", {expand_args}),')

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


def _output_to_path(name, value, dag):
    if isinstance(value, Input):
        return repr(f'inputs/{value.name}')
    if isinstance(value, Literal):
        return repr(value.value)
    if isinstance(value, Field):
        if isinstance(value.source, StepCall):
            prev = value.source
            spec = (prev.task or {}).get('outputs', {}).get(value.name, {})
            default_spec = spec.get('default')
            if default_spec:
                rendered = _interp_to_smk(default_spec)
                if rendered is not None:
                    return repr(rendered)
            if prev.map is not None:
                scatter_ports = set(prev.map.get('scatter', []))
                return _default_output_path(prev.id, value.name, True, scatter_ports)
            return repr(f'results/{prev.id}/{value.name}')
        if isinstance(value.source, Input):
            return repr(f'inputs/{value.source.name}/{value.name}')
    if isinstance(value, Record):
        parts = []
        for fname, fval in sorted(value.fields.items()):
            parts.append(_output_to_path(f'{name}_{fname}', fval, dag))
        return ', '.join(parts)
    raise ValueError(f'Unsupported output value for Snakemake: {type(value).__name__}')


def _mapped_by_step_to_smk(step, dag):
    map_info = step.map or {}
    group_key = map_info.get('group_by')
    source = map_info.get('source', {})
    input_schema = step.input_schema or {}
    task = step.task or {}
    body = task.get('body', '')
    rname = _rule_name(step.id)
    grp_name = f'group_{rname}'
    key_wildcard = f'{_san(rname)}_key'
    lines = []

    col_names = list(input_schema.keys())
    if group_key and group_key in col_names:
        col_names = [group_key] + [n for n in col_names if n != group_key]

    src_name = source_input_name(source)
    src_cols = table_columns(source)

    lines.append(f'# map_by: group by {group_key}, then apply {rname}')
    lines.append('')

    lines.append(f'rule {grp_name}:')
    lines.append('    input:')
    if src_name:
        for col in col_names:
            col_binding = src_cols.get(col, {})
            col_src_name = col_binding.get('name', col) if isinstance(col_binding, dict) else col
            sep = '/' if col_src_name else ''
            pattern = 'inputs/' + col_src_name + sep + '{' + col + '}'
            lines.append(f'        {_san(col)}={pattern!r},')
        lines.append('    params:')
        lines.append(f'        group_key="""{group_key}""",')
    else:
        for col in col_names:
            binding = step.bindings.get(col)
            path_expr = _binding_to_path(binding, col, True, {col})
            lines.append(f'        {_san(col)}={path_expr},')
    lines.append('    output:')
    lines.append(f'        groups=directory("results/{rname}/group/{{{key_wildcard}}}"),')
    lines.append('    run:')
    lines.append('        import json, os')
    lines.append(f'        keys = set()')
    lines.append(f'        # TODO: read column data, group by {group_key}')
    lines.append(f'        # For each key value, write results to directory')
    lines.append(f'        os.makedirs(str(wildcards.{key_wildcard}), exist_ok=True)')
    lines.append('')

    lines.append(f'rule {rname}:')
    lines.append('    input:')
    lines.append(f'        groups="results/{rname}/group/{{{key_wildcard}}}",')
    lines.append('    output:')
    for out_name in step.outputs:
        default_spec = task.get('outputs', {}).get(out_name, {}).get('default')
        if default_spec:
            rendered = _interp_to_smk(default_spec)
            if rendered is not None:
                path_expr = repr(rendered)
            else:
                path_expr = repr(f'results/{{{key_wildcard}}}/{out_name}')
        else:
            path_expr = repr(f'results/{{{key_wildcard}}}/{out_name}')
        lines.append(f'        {_san(out_name)}={path_expr},')
    lines.append('    shell:')
    if body.strip():
        interp_body = _interpolate_shell(body, step)
        lines.append(f'        "{interp_body}"')
    else:
        lines.append(f'        "tool {{input.groups}} > {{output}}"')
    lines.append('')

    return '\n'.join(lines)


def _subworkflow_to_smk(step, parent_id):
    definition = step.task or {}
    dag_data = definition.get('dag', {})
    if not dag_data:
        return ''
    wf_id = f'{parent_id}_{step.id}'
    return transpile_dag_dict(dag_data, workflow_id=wf_id, _top_level=False)
