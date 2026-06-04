import json

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import emit_name, interp_script, run_value, source_kind, step_name, workflow_name, word_interp, _flatten_output_names
from swl.types import to_nf_qualifier


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def transpile_dag_dict(data, workflow_id='main', _top_level=True):
    dag = DAG.from_dict(data)
    dag.validate()
    dag.outputs = _flatten_output_names(dag.outputs)
    _validate_supported(dag)

    processes = {}
    sub_workflows = []
    for step in dag.steps:
        if step.id not in processes:
            if step.type == 'workflow':
                wdl = _subworkflow_to_nf(step, workflow_id)
                processes[step.id] = _process_name(step.id)
                sub_workflows.append(wdl)
            else:
                processes[step.id] = _process_name(step.id)

    lines = []
    for step in dag.steps:
        if step.type != 'workflow':
            lines.append(_task_to_process(step))
            lines.append('')
    for sw in sub_workflows:
        lines.append(sw)
        lines.append('')
    lines.append(_dag_to_nf(dag, workflow_id, processes))
    return '\n'.join(lines)



def _process_name(step_id):
    return step_name(step_id, 'PROCESS', upper=True)


def _nf_emit_name(name):
    return emit_name(name)


def _task_to_process(step):
    task = step.task or {}
    body = task.get('body', '')
    pname = _process_name(step.id)
    lines = [f'process {pname} {{', '']

    directives = _emit_directives(step)
    if directives:
        for d in directives.split('\n'):
            lines.append(f'    {d.strip()}')
        lines.append('')

    inputs = task.get('inputs', {})
    has_map = getattr(step, 'map', None) is not None
    if inputs:
        lines.append('    input:')
        if has_map:
            input_names = list(inputs.keys())
            tuple_parts = []
            for in_name in input_names:
                spec = inputs.get(in_name, {})
                qual, _ = to_nf_qualifier(spec.get('type'))
                if qual == 'path':
                    tuple_parts.append(f'path({in_name})')
                else:
                    tuple_parts.append(f'val({in_name})')
            lines.append(f'    tuple {"(" + ", ".join(tuple_parts) + ")" if len(tuple_parts) > 1 else tuple_parts[0]}')
        else:
            for in_name, spec in inputs.items():
                qual, _ = to_nf_qualifier(spec.get('type'))
                if qual == 'path':
                    lines.append(f'    path {in_name}')
                else:
                    lines.append(f'    val {in_name}')
        lines.append('')

    outputs = task.get('outputs', {})
    if outputs:
        lines.append('    output:')
        for out_name, spec in outputs.items():
            qual, _ = to_nf_qualifier(spec.get('type'))
            default = spec.get('default')
            if default:
                path_expr = _interp_to_nf(default)
                if qual == 'path':
                    lines.append(f'    path "{path_expr}", emit: {_nf_emit_name(out_name)}')
                else:
                    lines.append(f'    val "{path_expr}", emit: {_nf_emit_name(out_name)}')
            else:
                qualifier = 'path' if qual == 'path' else 'val'
                glob_pattern = f'*.{out_name}' if qual == 'path' else out_name
                lines.append(f'    {qualifier} "{glob_pattern}", emit: {_nf_emit_name(out_name)}')
        lines.append('')

    if body.strip():
        interp_body = _interpolate_shell(body, step)
        lines.append('    script:')
        lines.append('    """')
        for line in interp_body.split('\n'):
            lines.append(f'    {line}' if line.strip() else '')
        lines.append('    """')
        lines.append('')

    lines.append('}')
    return '\n'.join(lines)


def _emit_directives(step):
    directives = []
    run = (step.task or {}).get('run', {})
    cpu = run_value(run, 'cpu')
    memory = run_value(run, 'memory')
    time = run_value(run, 'time')
    image = run_value(run, 'image')
    if cpu is not None:
        directives.append(f'    cpus {cpu}')
    if memory is not None:
        directives.append(f"    memory '{memory} MB'")
    if time is not None:
        directives.append(f"    time '{time}m'")
    if image is not None:
        directives.append(f"    container '{image}'")
    return '\n'.join(directives)


def _interp_to_nf(value):
    return word_interp(value, lambda text: text, lambda name: f"${{{name}}}", lambda text: f"${{{text}}}")


_NF_RUN_MAP = {'cpu': 'cpus', 'memory': 'memory', 'time': 'time'}


def _interpolate_shell(body, step):
    task = step.task or {}
    input_names = set(task.get('inputs', {}).keys())
    run = task.get('run', {})
    run_names = set()
    for rv in _NF_RUN_MAP:
        if run.get(rv, {}).get('value') is not None:
            run_names.add(rv)
    nf_run = {k: _NF_RUN_MAP[k] for k in run_names}

    def var_fn(name):
        if name in input_names:
            return '${' + name + '}'
        if name in nf_run:
            return '${task.' + nf_run[name] + '}'
        return None

    def expr_fn(text):
        resolved = text
        if 'memory' in resolved and 'memory' in nf_run:
            resolved = resolved.replace('memory', 'task.memory', 1)
        if 'cpu' in resolved and 'cpu' in nf_run:
            resolved = resolved.replace('cpu', 'task.cpus', 1)
        if 'time' in resolved and 'time' in nf_run:
            resolved = resolved.replace('time', 'task.time', 1)
        return '${' + resolved + '}'

    return interp_script(body, var_fn, expr_fn, joiner='')


def _channel_name(name):
    return f'{name}_ch'


def _input_channel(name, spec):
    swl_type = spec.type if hasattr(spec, 'type') else spec.get('type')
    if swl_type == 'file':
        return f'Channel.fromPath(params.{name}, checkIfExists: true)'
    if swl_type and swl_type.startswith('['):
        return f'Channel.fromPath(params.{name}, checkIfExists: true).toList()'
    return f'Channel.value(params.{name})'


def _dag_to_nf(dag, workflow_id, processes):
    lines = ['workflow {']

    channels = {}
    for name, spec in dag.inputs.items():
        ch_name = _channel_name(name)
        channels[name] = ch_name
        lines.append(f'    {ch_name} = {_input_channel(name, spec)}')

    if dag.inputs:
        lines.append('')

    for step in dag.steps:
        if getattr(step, 'map', None) is not None:
            map_info = step.map
            if map_info.get('group_by') is not None:
                lines.extend(_mapped_by_step_to_call(step, channels, processes))
            else:
                lines.extend(_mapped_step_to_call(step, channels, processes))
            lines.append('')
            continue

        if step.type == 'workflow':
            pname = _wf_name(f'{workflow_id}_{step.id}')
        else:
            pname = processes.get(step.id, _process_name(step.id))

        inputs = []
        for input_name in (step.task or {}).get('inputs', {}).keys():
            binding = step.bindings.get(input_name)
            if binding is None:
                ch = channels.get(input_name, _channel_name(input_name))
                inputs.append(ch)
            else:
                ch = _binding_to_channel(binding, channels, step)
                inputs.append(ch)

        if len(inputs) == 1:
            lines.append(f'    {pname}({inputs[0]})')
        elif len(inputs) > 1:
            ch_var = f'{step.id}_ch'
            join_expr = None
            for ch in inputs:
                join_expr = ch if join_expr is None else f'{join_expr}.join({ch})'
            lines.append(f'    {ch_var} = {join_expr}')
            lines.append(f'    {pname}({ch_var})')
        else:
            lines.append(f'    {pname}()')

        for out_name in step.outputs:
            channels[f'{step.id}.{out_name}'] = f'{pname}.out.{_nf_emit_name(out_name)}'

    if dag.outputs:
        lines.append('')
        for name, output in dag.outputs.items():
            binding = output.value if isinstance(output, OutputSpec) else output
            ch_expr = _binding_to_channel(binding, channels, None)
            lines.append(f'    {name}_out = {ch_expr}')
            lines.append(f'    emit: {name}_out')

    lines.append('}')
    return '\n'.join(lines)


def _binding_to_channel(binding, channels, current_step):
    if isinstance(binding, Input):
        return channels[binding.name]

    if isinstance(binding, Literal):
        val = json.dumps(binding.value)
        return f'Channel.value({val})'

    if isinstance(binding, Field):
        source_ch = _binding_to_channel(binding.source, channels, current_step)
        if isinstance(binding.source, StepCall):
            pname = _process_name(binding.source.id)
            ch = f'{pname}.out.{_nf_emit_name(binding.name)}'
            if getattr(binding.source, 'map', None) is not None:
                ch += '.toList()'
            return ch
        return f'{source_ch}.map{{ it.{binding.name} }}'

    if isinstance(binding, Merge):
        left = _binding_to_channel(binding.left, channels, current_step)
        right = _binding_to_channel(binding.right, channels, current_step)
        return f'{left}.join({right})'

    if isinstance(binding, Record):
        raise ValueError(
            f'Record binding with fields {list(binding.fields.keys())} must be flattened '
            'before Nextflow transpilation. This typically means a non-saturating record '
            'survived to DAG output. Use explicit field bindings instead.'
        )

    if isinstance(binding, StepCall):
        pname = _process_name(binding.id)
        return f'{pname}.out'

    raise ValueError(f'Unsupported binding for Nextflow: {type(binding).__name__}')


def _mapped_step_to_call(step, channels, processes):
    pname = _process_name(step.id)
    lines = []

    task_inputs = (step.task or {}).get('inputs', {})
    join_expr = None
    for input_name in task_inputs.keys():
        binding = step.bindings.get(input_name)
        if binding is None:
            ch = channels.get(input_name)
            if ch is None:
                if input_name in (getattr(step, 'input_schema', None) or {}):
                    ch = _channel_name(input_name)
                    channels[input_name] = ch
            if ch is None:
                continue
        else:
            ch = _binding_to_channel(binding, channels, step)
        join_expr = ch if join_expr is None else f'{join_expr}.join({ch})'

    ch_var = f'{step.id}_ch'
    if join_expr is not None:
        lines.append(f'    {ch_var} = {join_expr}')
        lines.append(f'    {pname}({ch_var})')
    else:
        lines.append(f'    {pname}()')

    for out_name in step.outputs:
        channels[f'{step.id}.{out_name}'] = f'{pname}.out.{_nf_emit_name(out_name)}'

    return lines


def _mapped_by_step_to_call(step, channels, processes):
    map_info = step.map or {}
    source = map_info.get('source', {})
    group_key = map_info.get('group_by')
    pname = _process_name(step.id)

    source_type = source_kind(source)
    src_ch = None
    if source_type == 'input':
        src_ch = channels.get(source.get('name')) or _channel_name(source['name'])

    if source_type == 'input' and src_ch is not None:
        group_key_nf = group_key.replace('-', '_')
        ch_var = f'{step.id}_records'
        lines = [f'    {ch_var} = {src_ch}.flatten()']
        grouped_var = f'{step.id}_groups'
        lines.append(f'    {grouped_var} = {ch_var}')
        lines.append(f'        .map{{ rec -> tuple(groupKey(rec.{group_key_nf}), rec) }}')
        lines.append(f'        .groupTuple()')
    else:
        ch_var = f'{step.id}_vals'
        if source_type == 'table':
            col_binding = source.get('columns', {}).get(group_key, {})
            col_src = col_binding.get('name', group_key)
            src_ch = channels.get(col_src) or _channel_name(col_src)
            lines = [f'    {ch_var} = {src_ch}']
        else:
            lines = [f'    {ch_var} = {src_ch or "Channel.empty()"}']
        grouped_var = f'{step.id}_groups'
        lines.append(f'    {grouped_var} = {ch_var}')
        lines.append(f'        .map{{ val -> tuple(groupKey(val), val) }}')
        lines.append(f'        .groupTuple()')

    lines.append('')
    lines.append(f'    {pname}({grouped_var})')

    for out_name in step.outputs:
        channels[f'{step.id}.{out_name}'] = f'{pname}.out.{_nf_emit_name(out_name)}'

    return lines


def _subworkflow_to_nf(step, parent_id):
    definition = step.task or {}
    dag_data = definition.get('dag', {})
    if not dag_data:
        return ''
    wf_id = f'{parent_id}_{step.id}'
    return transpile_dag_dict(dag_data, workflow_id=wf_id, _top_level=False)


def _wf_name(workflow_id):
    return workflow_name(workflow_id, 'main', upper=True)


def _validate_supported(dag):
    for name, output in dag.outputs.items():
        binding = output.value if isinstance(output, OutputSpec) else output
        if isinstance(binding, Literal):
            raise ValueError(f'Nextflow does not support literal workflow outputs: {name}')
