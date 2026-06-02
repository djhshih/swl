import json

from swl.ir.dag import DAG, Field, Input, Literal, Merge, Record, StepCall


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def transpile_dag_dict(data, workflow_id='main', _top_level=True):
    dag = DAG.from_dict(data)
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


def _input_qualifier(swl_type):
    return {
        'file': ('path', None),
        'str': ('val', 'string'),
        'int': ('val', 'integer'),
        'float': ('val', 'float'),
    }.get(swl_type, ('val', 'string'))


def _process_name(step_id):
    name = step_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'PROCESS'
    if name[0].isdigit():
        name = '_' + name
    return name.upper()


def _nf_emit_name(name):
    return name.replace('-', '_')


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
                qual, _ = _input_qualifier(spec.get('type'))
                if qual == 'path':
                    tuple_parts.append(f'path({in_name})')
                else:
                    tuple_parts.append(f'val({in_name})')
            lines.append(f'    tuple {"(" + ", ".join(tuple_parts) + ")" if len(tuple_parts) > 1 else tuple_parts[0]}')
        else:
            for in_name, spec in inputs.items():
                qual, _ = _input_qualifier(spec.get('type'))
                if qual == 'path':
                    lines.append(f'    path {in_name}')
                else:
                    lines.append(f'    val {in_name}')
        lines.append('')

    outputs = task.get('outputs', {})
    if outputs:
        lines.append('    output:')
        for out_name, spec in outputs.items():
            qual, _ = _input_qualifier(spec.get('type'))
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
        lines.append('    script:')
        lines.append('    """')
        for line in body.split('\n'):
            lines.append(f'    {line}' if line.strip() else '')
        lines.append('    """')
        lines.append('')

    lines.append('}')
    return '\n'.join(lines)


def _emit_directives(step):
    directives = []
    task = step.task or {}
    run = task.get('run', {})

    for name, spec in run.items():
        if not isinstance(spec, dict):
            continue
        value = spec.get('value')
        if value is None:
            continue

        if name == 'cpu':
            directives.append(f'    cpus {value}')
        elif name == 'memory':
            directives.append(f"    memory '{value} MB'")
        elif name == 'time':
            directives.append(f"    time '{value}m'")
        elif name == 'image':
            directives.append(f"    container '{value}'")

    return '\n'.join(directives)


def _interp_to_nf(value):
    if value is None:
        return None
    if value.get('kind') == 'word':
        parts = value.get('parts', [])
        result = ''
        for part in parts:
            if part.get('kind') == 'literal':
                result += part['text']
            elif part.get('kind') == 'var':
                result += f"${{{part['name']}}}"
            elif part.get('kind') == 'expr':
                result += f"${{{part['text']}}}"
        return result
    return None


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
        for name, binding in dag.outputs.items():
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
            'Record bindings must be flattened before Nextflow transpilation'
        )

    if isinstance(binding, dict):
        mapped_steps = _collect_mapped_steps(current_step)
        return _dict_binding_to_channel(binding, channels, current_step, mapped_steps)

    if isinstance(binding, StepCall):
        pname = _process_name(binding.id)
        return f'{pname}.out'

    raise ValueError(f'Unsupported binding for Nextflow: {type(binding).__name__}')


def _dict_binding_to_channel(binding, channels, current_step, mapped_steps=None):
    source = binding.get('source')
    if source == 'input':
        return channels[binding['name']]
    if source == 'literal':
        val = json.dumps(binding.get('value'))
        return f'Channel.value({val})'
    if source == 'field':
        inner = _dict_binding_to_channel(binding['value'], channels, current_step, mapped_steps)
        return f'{inner}.map{{ it.{binding["field"]} }}'
    if 'step' in binding and 'output' in binding:
        pname = _process_name(binding['step'])
        ch = f'{pname}.out.{_nf_emit_name(binding["output"])}'
        if mapped_steps and binding['step'] in mapped_steps:
            ch += '.toList()'
        return ch
    if source == 'record':
        raise ValueError(
            'Record bindings must be flattened before Nextflow transpilation'
        )
    if source == 'merge':
        raise ValueError(
            'Merge bindings must be flattened before Nextflow transpilation'
        )
    raise ValueError(f'Unsupported binding source for Nextflow: {source!r}')


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

    source_type = source.get('source') if isinstance(source, dict) else None
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
    name = workflow_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'main'
    if name[0].isdigit():
        name = '_' + name
    return name.upper()


def _collect_mapped_steps(step):
    if step is None:
        return set()
    mapped = set()
    for name, binding in step.bindings.items():
        if isinstance(binding, Field) and isinstance(binding.source, StepCall):
            if getattr(binding.source, 'map', None) is not None:
                mapped.add(binding.source.id)
    return mapped


def _validate_supported(dag):
    for step in dag.steps:
        for name, binding in step.bindings.items():
            _validate_binding(binding, step)

    for name, binding in dag.outputs.items():
        _validate_output_binding(binding, name)


def _validate_binding(value, step):
    if isinstance(value, Merge):
        raise ValueError(
            f'Nextflow does not support merge bindings '
            f'(step {step.id}). Flatten merges before transpiling.'
        )
    if isinstance(value, dict) and value.get('source') == 'merge':
        raise ValueError(
            f'Nextflow does not support merge bindings '
            f'(step {step.id}). Flatten merges before transpiling.'
        )


def _validate_output_binding(value, name):
    if isinstance(value, Literal):
        raise ValueError(f'Nextflow does not support literal workflow outputs: {name}')
    if isinstance(value, Merge):
        raise ValueError(f'Nextflow does not support merge workflow outputs: {name}')
