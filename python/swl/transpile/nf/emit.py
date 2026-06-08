import json

from swl.dag.node import DAG, Field, Input, Literal, OutputSpec, Record, StepCall
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
    _collect_processes(dag, processes)

    lines = []
    _emit_processes(dag, lines, processes)
    lines.append(_dag_to_nf(dag, workflow_id, processes))
    return '\n'.join(lines)


def _collect_processes(dag, processes):
    for step in dag.steps:
        if step.type != 'workflow':
            processes[step.id] = _process_name(step.id)
        else:
            inner = step.task or {}
            inner_dag = inner.get('dag', {})
            if inner_dag:
                _collect_processes(DAG.from_dict(inner_dag), processes)


def _emit_processes(dag, lines, processes):
    for step in dag.steps:
        if step.type != 'workflow':
            lines.append(_task_to_process(step))
            lines.append('')
        else:
            inner = step.task or {}
            inner_dag = inner.get('dag', {})
            if inner_dag:
                _emit_processes(DAG.from_dict(inner_dag), lines, processes)



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


def _inline_dag_steps(inner_dag, channels, lines, processes):
    for step in inner_dag.steps:
        if getattr(step, 'map', None) is not None:
            map_info = step.map
            if map_info.get('group_by') is not None:
                lines.extend(_mapped_by_step_to_call(step, channels, processes))
            else:
                lines.extend(_mapped_step_to_call(step, channels, processes))
            lines.append('')
            continue

        if step.type == 'workflow':
            _inline_workflow_step(step, channels, lines, processes)
            continue

        pname = _process_name(step.id)

        inputs = []
        for input_name in (step.task or {}).get('inputs', {}).keys():
            binding = step.bindings.get(input_name)
            if binding is None:
                ch = channels.get(input_name, _channel_name(input_name))
                inputs.append(ch)
            else:
                ch = _binding_to_channel(binding, channels, step)
                inputs.append(ch)

        if len(inputs) >= 1:
            lines.append(f'    {pname}({", ".join(inputs)})')
        else:
            lines.append(f'    {pname}()')

        for out_name in step.outputs:
            channels[f'{step.id}.{out_name}'] = f'{pname}.out.{_nf_emit_name(out_name)}'


def _inline_workflow_step(step, channels, lines, processes):
    definition = step.task or {}
    dag_data = definition.get('dag', {})
    if not dag_data:
        return
    inner_dag = DAG.from_dict(dag_data)
    _setup_inner_channels(step, inner_dag, channels)
    _inline_dag_steps(inner_dag, channels, lines, processes)
    for out_name, output in inner_dag.outputs.items():
        binding = output.value if isinstance(output, OutputSpec) else output
        ch = _binding_to_channel(binding, channels, None)
        channels[f'{step.id}.{out_name}'] = ch


def _setup_inner_channels(step, inner_dag, channels):
    for name in inner_dag.inputs:
        if name not in channels:
            binding = step.bindings.get(name)
            if binding is not None:
                ch = _binding_to_channel(binding, channels, step)
                channels[name] = ch


def _dag_to_nf(dag, workflow_id, processes):
    wf_header = f'workflow {_wf_name(workflow_id)} {{' if workflow_id != 'main' else 'workflow {'
    lines = [wf_header]

    input_sourced_sources = set()
    for step in dag.steps:
        m = getattr(step, 'map', None) or {}
        src = m.get('source', {})
        if isinstance(src, dict) and src.get('source') == 'input':
            input_sourced_sources.add(src['name'])

    channels = {}
    for name, spec in dag.inputs.items():
        ch_name = _channel_name(name)
        channels[name] = ch_name
        if name in input_sourced_sources:
            lines.append(f'    {ch_name} = Channel.fromList(params.{name})')
        else:
            lines.append(f'    {ch_name} = {_input_channel(name, spec)}')

    if dag.inputs:
        lines.append('')

    _inline_dag_steps(dag, channels, lines, processes)

    if dag.outputs:
        lines.append('')
        for name, output in dag.outputs.items():
            binding = output.value if isinstance(output, OutputSpec) else output
            ch_expr = _binding_to_channel(binding, channels, None)
            lines.append(f'    {name} = {ch_expr}')
            lines.append(f'    emit: {name}')

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
            step_key = f'{binding.source.id}.{binding.name}'
            if step_key in channels:
                return channels[step_key]
            pname = _process_name(binding.source.id)
            ch = f'{pname}.out.{_nf_emit_name(binding.name)}'
            if getattr(binding.source, 'map', None) is not None:
                ch += '.toList()'
            return ch
        return f'{source_ch}.map{{ it.{binding.name} }}'

    if isinstance(binding, Record):
        raise ValueError(
            f'Record binding with fields {list(binding.fields.keys())} must be flattened '
            'before Nextflow transpilation. This typically means a non-saturating record '
            'survived to DAG output. Use explicit field bindings instead.'
        )

    if isinstance(binding, StepCall):
        pname = _process_name(binding.id)
        if binding.id in channels:
            return channels[binding.id]
        return f'{pname}.out'

    raise ValueError(f'Unsupported binding for Nextflow: {type(binding).__name__}')


def _mapped_workflow_step_to_call(step, channels, processes):
    lines = []
    m = step.map
    schema = getattr(step, 'input_schema', None) or {}
    src = m.get('source', {})
    src_type = src.get('source')

    definition = step.task or {}
    dag_data = definition.get('dag', {})
    inner_dag = DAG.from_dict(dag_data) if dag_data else None

    if inner_dag is None:
        return lines

    field_names = list(inner_dag.inputs.keys())

    if src_type == 'input':
        src_name = src['name']
        base_ch = channels.get(src_name, _channel_name(src_name))
        for field_name in field_names:
            wdl_type = schema.get(field_name, 'str')
            qual, _ = to_nf_qualifier(wdl_type)
            if qual == 'path':
                lines.append(f'    {field_name} = {base_ch}.map{{ x -> file(x.{field_name}) }}')
            else:
                lines.append(f'    {field_name} = {base_ch}.map{{ x -> x.{field_name} }}')
        lines.append('')
        inner_channels = {name: name for name in field_names}
    elif src_type == 'table':
        columns = src.get('columns', {})
        inner_channels = {}
        for field_name in field_names:
            col = columns.get(field_name, {})
            if col.get('source') == 'input':
                outer_name = col['name']
                inner_channels[field_name] = channels.get(outer_name, _channel_name(outer_name))
            else:
                inner_channels[field_name] = _channel_name(field_name)
    else:
        return lines

    _inline_dag_steps(inner_dag, inner_channels, lines, processes)

    for out_name in step.outputs:
        inner_output_binding = inner_dag.outputs.get(out_name)
        if inner_output_binding is not None:
            binding = inner_output_binding.value if isinstance(inner_output_binding, OutputSpec) else inner_output_binding
            ch = _binding_to_channel(binding, inner_channels, None)
        else:
            ch = f'{_process_name(step.id)}.out.{_nf_emit_name(out_name)}'
        channels[f'{step.id}.{out_name}'] = f'{ch}.toList()'

    return lines


def _mapped_step_to_call(step, channels, processes):
    pname = _process_name(step.id)
    lines = []

    m = getattr(step, 'map', None) or {}
    src = m.get('source', {})

    if step.type == 'workflow' and isinstance(src, dict) and src.get('source') in ('input', 'table') and not m.get('group_by'):
        return _mapped_workflow_step_to_call(step, channels, processes)

    if isinstance(src, dict) and src.get('source') == 'input' and not m.get('group_by'):
        src_name = src['name']
        schema = getattr(step, 'input_schema', None) or {}
        step_inputs = getattr(step, 'inputs', None) or {}
        base_ch = channels.get(src_name, _channel_name(src_name))
        ch_var = f'{step.id}_ch'
        order = list(step_inputs.keys()) or sorted(schema)
        tuple_parts = []
        for col_name in order:
            wdl_type = schema.get(col_name, 'str')
            qual, _ = to_nf_qualifier(wdl_type)
            if qual == 'path':
                tuple_parts.append(f'file(x.{col_name})')
            else:
                tuple_parts.append(f'x.{col_name}')
        lines.append(f'    {ch_var} = {base_ch}')
        lines.append(f'        .map{{ x -> tuple(')
        for i, part in enumerate(tuple_parts):
            comma = ',' if i < len(tuple_parts) - 1 else ''
            lines.append(f'            {part}{comma}')
        lines.append(f'        ) }}')
        lines.append(f'    {pname}({ch_var})')
    else:
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
        schema = getattr(step, 'input_schema', None) or {}
        step_inputs = getattr(step, 'inputs', None) or {}
        group_key_nf = group_key.replace('-', '_')
        base_ch = channels.get(source.get('name'), _channel_name(source['name']))
        ch_var = f'{step.id}_ch'
        order = list(step_inputs.keys()) or sorted(schema)
        tuple_parts = []
        if group_key in order:
            wdl_type = schema.get(group_key, 'str')
            qual, _ = to_nf_qualifier(wdl_type)
            if qual == 'path':
                tuple_parts.append(f'file(x.{group_key_nf})')
            else:
                tuple_parts.append(f'x.{group_key_nf}')
        for col_name in order:
            if col_name == group_key:
                continue
            wdl_type = schema.get(col_name, 'str')
            qual, _ = to_nf_qualifier(wdl_type)
            if qual == 'path':
                tuple_parts.append(f'file(x.{col_name})')
            else:
                tuple_parts.append(f'x.{col_name}')
        lines = [f'    {ch_var} = {base_ch}']
        lines.append(f'        .map{{ x -> tuple(')
        for i, part in enumerate(tuple_parts):
            comma = ',' if i < len(tuple_parts) - 1 else ''
            lines.append(f'            {part}{comma}')
        lines.append(f'        ) }}')
        lines.append(f'        .groupTuple()')
        grouped_var = ch_var
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


def _wf_name(workflow_id):
    return workflow_name(workflow_id, 'main', upper=True)


def _validate_supported(dag):
    for name, output in dag.outputs.items():
        binding = output.value if isinstance(output, OutputSpec) else output
        if isinstance(binding, Literal):
            raise ValueError(f'Nextflow does not support literal workflow outputs: {name}')
