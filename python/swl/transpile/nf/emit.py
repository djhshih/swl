import json

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import emit_name, flatten_dag_outputs, format_resource_directives, interp_script, source_kind, step_name, validate_dag_for_transpile, workflow_name, word_interp
from swl.types import to_nf_qualifier


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def transpile_dag_dict(data, workflow_id='main', _top_level=True):
    dag = DAG.from_dict(data)
    dag.validate()
    flatten_dag_outputs(dag)
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
            inner_dag = step.task_def.get('dag', {})
            if inner_dag:
                _collect_processes(DAG.from_dict(inner_dag), processes)


def _emit_processes(dag, lines, processes):
    for step in dag.steps:
        if step.type != 'workflow':
            lines.append(_task_to_process(step))
            lines.append('')
        else:
            inner_dag = step.task_def.get('dag', {})
            if inner_dag:
                _emit_processes(DAG.from_dict(inner_dag), lines, processes)



def _process_name(step_id):
    return step_name(step_id, 'PROCESS', upper=True)


def _nf_emit_name(name):
    return emit_name(name)


def _task_to_process(step):
    task = step.task_def
    body = task.get('body', '')
    pname = _process_name(step.id)
    lines = [f'process {pname} {{', '']

    directives = _emit_directives(step)
    if directives:
        for d in directives.split('\n'):
            lines.append(f'    {d.strip()}')
        lines.append('')

    inputs = task.get('inputs', {})
    outputs = task.get('outputs', {})
    has_map = step.is_mapped

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
    run = step.task_def.get('run', {})
    directives = format_resource_directives(run, {
        'cpu': lambda v: f'    cpus {v}',
        'memory': lambda v: f"    memory '{v} MB'",
        'time': lambda v: f"    time '{v}m'",
        'image': lambda v: f"    container '{v}'",
    })
    return '\n'.join(directives)


def _interp_to_nf(value):
    return word_interp(value, lambda text: text, lambda name: f"${{{name}}}", lambda text: f"${{{text}}}")


_NF_RUN_MAP = {'cpu': 'cpus', 'memory': 'memory', 'time': 'time'}

def _interpolate_shell(body, step):
    if '$' not in body:
        return body
    task = step.task_def
    input_names = set(task.get('inputs', {}).keys())
    run = task.get('run', {})
    run_names = set()
    for rv in ('cpu', 'memory', 'time'):
        if run.get(rv, {}).get('value') is not None:
            run_names.add(rv)
    _NF_RUN_MAP = {'cpu': 'cpus', 'memory': 'memory', 'time': 'time'}

    def var_fn(name):
        if name in input_names:
            return f'${{{name}}}'
        if name in run_names:
            return f'${{task.{_NF_RUN_MAP.get(name, name)}}}'
        return None

    def expr_fn(text):
        resolved = text
        for swl_name, nf_name in _NF_RUN_MAP.items():
            if swl_name in run_names and swl_name in resolved:
                resolved = resolved.replace(swl_name, f'task.{nf_name}', 1)
        return '${' + resolved + '}'

    return interp_script(body, var_fn, expr_fn)


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
        if step.is_mapped:
            if step.has_group_by:
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
        for input_name in step.task_def.get('inputs', {}).keys():
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
    dag_data = step.task_def.get('dag', {})
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
    is_entrypoint = workflow_id == 'main'
    wf_header = f'workflow {_wf_name(workflow_id)} {{' if workflow_id != 'main' else 'workflow {'
    lines = [wf_header]

    input_sourced_sources = set()
    table_sourced_inputs = set()
    for step in dag.steps:
        src = step.map_source
        if src.get('source') == 'input':
            input_sourced_sources.add(src['name'])
        if src.get('source') == 'table':
            columns = src.get('columns', {})
            for col_binding in columns.values():
                if col_binding.get('source') == 'input':
                    table_sourced_inputs.add(col_binding['name'])

    channels = {}
    all_table = bool(dag.inputs) and set(dag.inputs.keys()) == table_sourced_inputs

    if all_table:
        lines.append('    xs_ch = Channel.fromList(params.xs)')
        for name in dag.inputs:
            ch_name = _channel_name(name)
            channels[name] = ch_name
            lines.append(f'    {ch_name} = xs_ch.map{{ it.{name} }}')
    else:
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
        if not is_entrypoint:
            lines.append('')
            lines.append('    emit:')
            for name in dag.outputs:
                lines.append(f'    {name}')

    lines.append('}')
    return '\n'.join(lines)


def _binding_to_channel(binding, channels, current_step):
    if isinstance(binding, Input):
        return channels[binding.name]
    if isinstance(binding, Literal):
        val = json.dumps(binding.value)
        return f'Channel.value({val})'
    if isinstance(binding, Field):
        if isinstance(binding.source, StepCall):
            step_key = f'{binding.source.id}.{binding.name}'
            if step_key in channels:
                return channels[step_key]
            pname = _process_name(binding.source.id)
            return f'{pname}.out.{_nf_emit_name(binding.name)}'
        return channels.get(binding.name, _channel_name(binding.name))
    if isinstance(binding, Record):
        raise ValueError(
            f'Record binding with fields {list(binding.fields.keys())} must be flattened '
            'before Nextflow transpilation.'
        )
    if isinstance(binding, Merge):
        raise ValueError(
            f'Merge bindings must be flattened before Nextflow transpilation.'
        )
    if isinstance(binding, StepCall):
        if binding.id in channels:
            return channels[binding.id]
        pname = _process_name(binding.id)
        return f'{pname}.out'
    raise ValueError(f'Unsupported binding for Nextflow: {type(binding).__name__}')


def _mapped_workflow_step_to_call(step, channels, processes):
    lines = []
    src = step.map_source
    src_type = src.get('source')
    schema = step.input_schema_def

    dag_data = step.task_def.get('dag', {})
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

    src = step.map_source

    if step.type == 'workflow' and src.get('source') in ('input', 'table') and not step.map_info.get('group_by'):
        return _mapped_workflow_step_to_call(step, channels, processes)

    if src.get('source') == 'input' and not step.map_info.get('group_by'):
        src_name = src['name']
        schema = step.input_schema_def
        base_ch = channels.get(src_name, _channel_name(src_name))
        ch_var = f'{step.id}_ch'
        order = sorted(schema)
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
        task_inputs = step.task_def.get('inputs', {})
        join_expr = None
        for input_name in task_inputs.keys():
            binding = step.bindings.get(input_name)
            if binding is None:
                ch = channels.get(input_name)
                if ch is None:
                    if input_name in step.input_schema_def:
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
    source = step.map_source
    group_key = step.map_info.get('group_by')

    if step.type == 'workflow' and source.get('source') in ('input', 'table'):
        if step.has_group_by and source.get('source') == 'table':
            return _mapped_by_workflow_step(step, channels, processes, source, group_key)
        return _mapped_workflow_step_to_call(step, channels, processes)

    pname = _process_name(step.id)

    source_type = source_kind(source)
    src_ch = None
    if source_type == 'input':
        src_ch = channels.get(source.get('name')) or _channel_name(source['name'])

    if source_type == 'input' and src_ch is not None:
        schema = step.input_schema_def
        group_key_nf = group_key.replace('-', '_')
        base_ch = channels.get(source.get('name'), _channel_name(source['name']))
        ch_var = f'{step.id}_ch'
        order = sorted(schema)
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


def _mapped_by_workflow_step(step, channels, processes, source, group_key):
    source_type = source_kind(source)
    lines = []

    dag_data = step.task_def.get('dag', {})
    inner_dag = DAG.from_dict(dag_data) if dag_data else None
    if inner_dag is None:
        return lines

    field_names = list(inner_dag.inputs.keys())

    if source_type == 'input':
        src_name = source['name']
        base_ch = channels.get(src_name, _channel_name(src_name))
        schema = step.input_schema_def
        group_key_nf = group_key.replace('-', '_')
        order = sorted(schema)
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
        ch_var = f'{step.id}_ch'
        lines.append(f'    {ch_var} = {base_ch}')
        lines.append(f'        .map{{ x -> tuple(')
        for i, part in enumerate(tuple_parts):
            comma = ',' if i < len(tuple_parts) - 1 else ''
            lines.append(f'            {part}{comma}')
        lines.append(f'        ) }}')
        lines.append(f'        .groupTuple()')

        inner_channels = {}
        for i, field_name in enumerate(field_names):
            inner_channels[field_name] = f'{ch_var}.map{{ it[{(i + 1) if field_name == group_key else i}] }}'
    elif source_type == 'table':
        columns = source.get('columns', {})
        ch_var = f'{step.id}_vals'
        col_binding = columns.get(group_key, {})
        col_src = col_binding.get('name', group_key)
        src_ch = channels.get(col_src) or _channel_name(col_src)
        lines.append(f'    {ch_var} = {src_ch}')
        grouped_var = f'{step.id}_groups'
        lines.append(f'    {grouped_var} = {ch_var}')
        lines.append(f'        .map{{ val -> tuple(groupKey(val), val) }}')
        lines.append(f'        .groupTuple()')

        inner_channels = {}
        for field_name in field_names:
            col = columns.get(field_name, {})
            if col.get('source') == 'input' and field_name == group_key:
                inner_channels[field_name] = f'{grouped_var}.map{{ it[1] }}'
            elif col.get('source') == 'input':
                inner_channels[field_name] = channels.get(col['name'], _channel_name(col['name']))
            else:
                inner_channels[field_name] = _channel_name(field_name)
    else:
        return lines

    if lines:
        lines.append('')
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


def _wf_name(workflow_id):
    return workflow_name(workflow_id, 'main', upper=True)


def _validate_supported(dag):
    validate_dag_for_transpile(dag, 'Nextflow')
