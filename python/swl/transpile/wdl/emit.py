import json

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import (
    column_input_name,
    interp_script,
    source_input_name,
    source_kind,
    step_name,
    table_columns,
    workflow_name,
    word_interp,
    _flatten_output_names,
)
from swl.types import to_wdl_type

_RUN_VAR_TYPE = {'cpu': 'Int', 'memory': 'Int', 'time': 'Int'}


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def _conflicting_outputs(step):
    task = step.task or {}
    inputs = task.get('inputs', {})
    outputs = task.get('outputs', {})
    renames = {}
    for name in outputs:
        if name in inputs:
            renames[name] = name + '_out'
    return renames


def _conflicting_outputs_from_data(step_data):
    inputs = step_data.get('inputs', {})
    outputs = step_data.get('outputs', {})
    renames = {}
    for name in outputs:
        if name in inputs:
            renames[name] = name + '_out'
    return renames


def _build_output_rename_map(dag):
    rename_map = {}
    for step in dag.steps:
        if step.type != 'workflow':
            renames = _conflicting_outputs(step)
            if renames:
                rename_map[step.id] = renames
    return rename_map


def transpile_dag_dict(data, workflow_id='main', _top_level=True):
    dag = DAG.from_dict(data)
    dag.validate()
    dag.outputs = _flatten_output_names(dag.outputs)
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if isinstance(value, Literal):
            raise ValueError(f'WDL does not support literal workflow outputs: {name}')

    structs = _collect_structs(dag)
    output_renames = _build_output_rename_map(dag)

    input_struct_map = {}
    for step in dag.steps:
        m = getattr(step, 'map', None) or {}
        src = m.get('source', {})
        if isinstance(src, dict) and src.get('source') == 'input' and not m.get('group_by'):
            src_name = src['name']
            if src_name not in input_struct_map:
                schema = step.input_schema or {}
                input_struct_map[src_name] = (
                    _emit_input_struct(step.id, schema),
                    _input_struct_name(step.id),
                )

    tasks = {}
    inline_map = {}
    inline_output_map = {}
    for step in dag.steps:
        if step.id not in tasks:
            step_map = getattr(step, 'map', None) or {}
            if step.type == 'workflow':
                sub_dag_data = (step.task or {}).get('dag', {})
                if sub_dag_data:
                    sub_steps_data = sub_dag_data.get('steps', [])
                    sub_renames = {}
                    for sd in sub_steps_data:
                        sid = sd['id']
                        if sid not in tasks:
                            tasks[sid] = _task_name(sid)
                        sr = _conflicting_outputs_from_data(sd)
                        if sr:
                            sub_renames[sid] = sr
                            if sid not in output_renames:
                                output_renames[sid] = sr
                    out_map = {}
                    for out_name, out_spec in sub_dag_data.get('outputs', {}).items():
                        val = (out_spec or {}).get('value', {}) if isinstance(out_spec, dict) else {}
                        src_step = val.get('step')
                        src_output = val.get('output')
                        if src_step and src_output:
                            sr = sub_renames.get(src_step, {})
                            out_map[out_name] = (_call_alias(src_step), sr.get(src_output, src_output))
                    inline_map[step.id] = sub_steps_data
                    if out_map:
                        inline_output_map[step.id] = out_map
            else:
                tasks[step.id] = _task_name(step.id)

    lines = []
    if _top_level:
        has_map_by = any(
            getattr(s, 'map', None) and s.map.get('group_by')
            for s in dag.steps
        )
        lines.append('version 1.1\n' if has_map_by else 'version 1.0\n')
    for s in structs:
        lines.append(s)
        lines.append('')
    for struct_def, _ in input_struct_map.values():
        lines.append(struct_def)
        lines.append('')
    for step in dag.steps:
        if step.type != 'workflow' and step.id not in [s.id for s in dag.steps if s.type == 'workflow']:
            lines.append(_task_to_wdl(step, output_renames.get(step.id)))
            lines.append('')
    emitted_sub_ids = set()
    for sub_steps in inline_map.values():
        for sd in sub_steps:
            sid = sd['id']
            if sid not in emitted_sub_ids:
                lines.append(_task_to_wdl_from_data(sid, sd))
                lines.append('')
                emitted_sub_ids.add(sid)
    lines.append(_dag_to_wdl(dag, workflow_id, tasks, input_struct_map, output_renames, inline_map, inline_output_map))
    return '\n'.join(lines)


_WDL_RESERVED = {
    'call', 'workflow', 'task', 'scatter', 'if', 'then', 'else',
    'while', 'input', 'output', 'command', 'runtime', 'meta',
    'parameter_meta', 'import', 'struct', 'version', 'as', 'in',
    'true', 'false', 'object', 'null',
    'Array', 'Map', 'Pair', 'String', 'File', 'Int', 'Float', 'Boolean', 'Object',
}


def _sanitize_ident(name):
    if name in _WDL_RESERVED:
        return name + '_'
    return name


def _task_name(step_id):
    return _sanitize_ident(step_name(step_id, 'task'))


def _wf_name(workflow_id):
    return workflow_name(workflow_id, 'main')


def _call_alias(step_id):
    return _sanitize_ident(step_name(step_id, 'step'))


def _task_to_wdl(step, rename_map=None):
    if rename_map is None:
        rename_map = {}
    task = step.task or {}
    body = task.get('body', '')
    lines = [f'task {_task_name(step.id)} {{', '']

    inputs = task.get('inputs', {})
    run = task.get('run', {})
    run_var_names = set()
    for rv in ('cpu', 'memory', 'time'):
        if run.get(rv, {}).get('value') is not None:
            run_var_names.add(rv)

    array_input_names = set()
    if inputs or run_var_names:
        lines.append('    input {')
        for name, spec in inputs.items():
            t = to_wdl_type(spec.get('type'))
            if t.startswith('Array['):
                array_input_names.add(name)
            lines.append(f'        {t} {name}')
        for rv in sorted(run_var_names):
            t = _RUN_VAR_TYPE[rv]
            val = run[rv]['value']
            lines.append(f'        {t} {rv} = {val}')
        lines.append('    }')
        lines.append('')

    lines.append('    command <<<')
    for line in body.split('\n'):
        interp_line = _interpolate_bash_vars(line, set(inputs.keys()) | run_var_names, array_input_names)
        lines.append(f'        {interp_line}')
    lines.append('    >>>')
    lines.append('')

    outputs = task.get('outputs', {})
    if outputs:
        lines.append('    output {')
        for name, spec in outputs.items():
            emit_name = rename_map.get(name, name)
            t = to_wdl_type(spec.get('type'))
            default = spec.get('default')
            if default:
                path_expr = _interp_to_wdl(default)
                lines.append(f'        {t} {emit_name} = "{path_expr}"')
            else:
                lines.append(f'        {t} {emit_name} = "*.{name}"')
        lines.append('    }')
        lines.append('')

    req = _emit_requirements(step)
    if req:
        lines.append(req)
        lines.append('')

    lines.append('}')
    return '\n'.join(lines)


def _task_to_wdl_from_data(step_id, step_data):
    rename_map = _conflicting_outputs_from_data(step_data)
    body = step_data.get('script', '')
    lines = [f'task {_task_name(step_id)} {{', '']

    inputs = step_data.get('inputs', {})
    run = step_data.get('run', {})
    run_var_names = set()
    for rv in ('cpu', 'memory', 'time'):
        if run.get(rv, {}).get('value') is not None:
            run_var_names.add(rv)

    array_input_names = set()
    if inputs or run_var_names:
        lines.append('    input {')
        for name, spec in inputs.items():
            t = to_wdl_type(spec.get('type'))
            if t.startswith('Array['):
                array_input_names.add(name)
            lines.append(f'        {t} {name}')
        for rv in sorted(run_var_names):
            t = _RUN_VAR_TYPE[rv]
            val = run[rv]['value']
            lines.append(f'        {t} {rv} = {val}')
        lines.append('    }')
        lines.append('')

    lines.append('    command <<<')
    for line in body.split('\n'):
        interp_line = _interpolate_bash_vars(line, set(inputs.keys()) | run_var_names, array_input_names)
        lines.append(f'        {interp_line}')
    lines.append('    >>>')
    lines.append('')

    outputs = step_data.get('outputs', {})
    if outputs:
        lines.append('    output {')
        for name, spec in outputs.items():
            emit_name = rename_map.get(name, name)
            t = to_wdl_type(spec.get('type'))
            default = spec.get('default')
            if default:
                path_expr = _interp_to_wdl(default)
                lines.append(f'        {t} {emit_name} = "{path_expr}"')
            else:
                lines.append(f'        {t} {emit_name} = "*.{name}"')
        lines.append('    }')
        lines.append('')

    attrs = []
    for name, spec in run.items():
        if not isinstance(spec, dict):
            continue
        value = spec.get('value')
        if value is None:
            continue
        if name == 'cpu':
            attrs.append(f'        cpu: {value}')
        elif name == 'memory':
            attrs.append(f'        memory: "{_format_memory(value)}"')
        elif name == 'image':
            attrs.append(f'        container: "{value}"')
        elif name == 'time':
            attrs.append(f'        time_minutes: {value}')
    if attrs:
        lines.append('    requirements {')
        lines.extend(attrs)
        lines.append('    }')
        lines.append('')

    lines.append('}')
    return '\n'.join(lines)


def _interpolate_bash_vars(body, known_vars, array_vars=None):
    if array_vars is None:
        array_vars = set()

    def var_fn(name):
        if name in known_vars:
            if name in array_vars:
                return '~{sep=" " ' + name + '}'
            return '~{' + name + '}'
        return None

    def expr_fn(text):
        return '~{' + text + '}'

    return interp_script(body, var_fn, expr_fn, joiner='')


def _emit_requirements(step):
    attrs = []
    task = step.task or {}
    run = task.get('run', {})

    for name, spec in run.items():
        if not isinstance(spec, dict):
            continue
        value = spec.get('value')
        if value is None:
            continue

        if name == 'cpu':
            attrs.append(f'        cpu: {value}')
        elif name == 'memory':
            memory_val = _format_memory(value)
            attrs.append(f'        memory: "{memory_val}"')
        elif name == 'image':
            attrs.append(f'        container: "{value}"')
        elif name == 'time':
            attrs.append(f'        time_minutes: {value}')

    if not attrs:
        return ''
    return '    requirements {\n' + '\n'.join(attrs) + '\n    }'


def _format_memory(value):
    if isinstance(value, str):
        return value
    return f'{value} MB'


def _interp_to_wdl(value):
    return word_interp(value, lambda text: text, lambda name: f"~{{{name}}}", lambda text: f"~{{{text}}}")


def _dag_to_wdl(dag, workflow_id, tasks, input_struct_map=None, output_renames=None, inline_map=None, inline_output_map=None):
    if input_struct_map is None:
        input_struct_map = {}
    if output_renames is None:
        output_renames = {}
    if inline_map is None:
        inline_map = {}
    if inline_output_map is None:
        inline_output_map = {}
    lines = [f'workflow {_wf_name(workflow_id)} {{', '']

    decomposed_inputs = {}
    input_blacklist = set()
    array_inputs = set()
    for step in dag.steps:
        m = getattr(step, 'map', None) or {}
        src = m.get('source', {})
        if isinstance(src, dict) and src.get('source') == 'input':
            src_name = src['name']
            input_blacklist.add(src_name)
            if m.get('group_by'):
                schema = step.input_schema or {}
                group_key = m.get('group_by')
                col_names = [group_key] + [n for n in schema if n != group_key]
                decomposed_inputs[step.id] = [
                    (col, f'Array[{to_wdl_type(schema.get(col, "str"))}]')
                    for col in col_names
                ]
        elif isinstance(src, dict) and src.get('source') == 'table':
            for col_name, col_src in (src.get('columns') or {}).items():
                if isinstance(col_src, dict) and col_src.get('source') == 'input':
                    array_inputs.add(col_src['name'])

    if dag.inputs:
        lines.append('    input {')
        for name, spec in dag.inputs.items():
            if name in input_struct_map:
                _, struct_name = input_struct_map[name]
                lines.append(f'        Array[{struct_name}] {name}')
            elif name in input_blacklist:
                continue
            elif name in array_inputs:
                t = to_wdl_type(spec.type)
                lines.append(f'        Array[{t}] {name}')
            else:
                t = to_wdl_type(spec.type)
                lines.append(f'        {t} {name}')
        for step_id, cols in decomposed_inputs.items():
            for col_name, col_type in cols:
                lines.append(f'        {col_type} {col_name}')
        lines.append('    }')
        lines.append('')

    for step in dag.steps:
        if getattr(step, 'map', None) is not None:
            if step.map.get('group_by') is not None:
                lines.extend(_mapped_by_step_to_wdl(step, tasks, output_renames, inline_map=inline_map))
            else:
                lines.extend(_mapped_step_to_wdl(step, tasks, output_renames, inline_map=inline_map))
            continue

        if step.id in inline_map:
            sub_steps = inline_map[step.id]
            for sd in sub_steps:
                sid = sd['id']
                tname = tasks.get(sid, _task_name(sid))
                alias = _call_alias(sid)
                if alias == tname:
                    lines.append(f'    call {tname} {{')
                else:
                    lines.append(f'    call {tname} as {alias} {{')
                lines.append('        input:')
                sub_inputs = sd.get('inputs', {})
                sub_bindings = sd.get('bindings', {})
                for in_name in sub_inputs:
                    binding = sub_bindings.get(in_name)
                    if binding is None:
                        expr = in_name
                    else:
                        expr = _sub_step_binding_to_expr(binding, output_renames)
                    lines.append(f'            {in_name} = {expr},')
                lines.append('    }')
                lines.append('')
            continue

        tname = tasks.get(step.id, _task_name(step.id))
        alias = _call_alias(step.id)
        if type(step).__name__ == 'StepCall' and step.type == 'workflow':
            tname = _wf_name(f'{workflow_id}_{step.id}')

        if alias == tname:
            lines.append(f'    call {tname} {{')
        else:
            lines.append(f'    call {tname} as {alias} {{')
        lines.append('        input:')
        for in_name, binding in step.bindings.items():
            expr = _binding_to_wdl_expr(binding, step.id, dag, output_renames, inline_output_map)
            lines.append(f'            {in_name} = {expr},')
        lines.append('    }')
        lines.append('')

    if dag.outputs:
        lines.append('    output {')
        for name, output in dag.outputs.items():
            binding = output.value if isinstance(output, OutputSpec) else output
            expr = _binding_to_wdl_expr(binding, None, dag, output_renames, inline_output_map)
            typ = to_wdl_type(output.type, output.optional)
            lines.append(f'        {typ} {name} = {expr}')
        lines.append('    }')

    lines.append('}')
    return '\n'.join(lines)


def _binding_to_wdl_expr(binding, current_step_id, dag, output_renames=None, inline_output_map=None):
    if output_renames is None:
        output_renames = {}
    if inline_output_map is None:
        inline_output_map = {}

    if isinstance(binding, Input):
        return binding.name

    if isinstance(binding, Literal):
        return _literal_to_wdl(binding.value)

    if isinstance(binding, Field):
        source_expr = _binding_to_wdl_expr(binding.source, current_step_id, dag, output_renames, inline_output_map)
        field_name = binding.name
        if isinstance(binding.source, StepCall):
            step_id = binding.source.id
            if step_id in inline_output_map and field_name in inline_output_map[step_id]:
                alias, out_name = inline_output_map[step_id][field_name]
                return f'{alias}.{out_name}'
            step_renames = output_renames.get(step_id, {})
            field_name = step_renames.get(field_name, field_name)
        return f'{source_expr}.{field_name}'

    if isinstance(binding, Merge):
        raise ValueError(
            f'Merge bindings must be flattened before WDL transpilation'
        )

    if isinstance(binding, Record):
        struct_name = _struct_name_for(binding)
        fields = []
        for fname, fbinding in binding.fields.items():
            fexpr = _binding_to_wdl_expr(fbinding, current_step_id, dag, output_renames, inline_output_map)
            fields.append(f'{fname}: {fexpr}')
        return f'{struct_name} {{{", ".join(fields)}}}'

    if isinstance(binding, StepCall):
        step_id = binding.id
        if step_id in inline_output_map:
            return _call_alias(step_id)
        return _call_alias(binding.id)

    raise ValueError(f'Unsupported binding for WDL: {type(binding).__name__}')


def _literal_to_wdl(value):
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if value is None:
        raise ValueError('None/null literals cannot be represented in WDL')
    return str(value)


def _sub_step_binding_to_expr(binding, output_renames=None):
    if output_renames is None:
        output_renames = {}
    src = binding.get('source')
    if src == 'input':
        return binding['name']
    if src == 'step_output':
        step_id = binding['step']
        out_name = binding['output']
        step_renames = output_renames.get(step_id, {})
        renamed = step_renames.get(out_name, out_name)
        return f'{_call_alias(step_id)}.{renamed}'
    if src == 'literal':
        return json.dumps(binding['value'])
    return binding.get('name', '???')


def _mapped_step_to_wdl(step, tasks, output_renames=None, inline_map=None):
    if output_renames is None:
        output_renames = {}
    if inline_map is None:
        inline_map = {}
    lines = []
    s_var = f'{step.id}_i'

    map_info = step.map or {}
    source = map_info.get('source', {})
    table_columns = _get_table_columns(source)
    len_expr = _derive_length_expr(source, step.bindings, table_columns)
    input_name = source_input_name(source)
    lines.append(f'    scatter ({s_var} in range({len_expr})) {{')

    if step.id in inline_map:
        sub_steps = inline_map[step.id]
        for sd in sub_steps:
            sid = sd['id']
            tname = tasks.get(sid, _task_name(sid))
            alias = _call_alias(sid)
            call_line = f'        call {tname}'
            if alias != tname:
                call_line += f' as {alias}'
            lines.append(call_line + ' {')
            lines.append('            input:')
            sub_inputs = sd.get('inputs', {})
            sub_bindings = sd.get('bindings', {})
            for in_name in sub_inputs:
                binding = sub_bindings.get(in_name)
                if binding is None:
                    if input_name:
                        expr = f'{input_name}[{s_var}].{in_name}'
                    elif in_name in table_columns:
                        expr = f'{_column_input_name(source, in_name)}[{s_var}]'
                    else:
                        expr = in_name
                elif binding.get('source') == 'step_output':
                    step_id = binding['step']
                    out_name = binding['output']
                    step_renames = output_renames.get(step_id, {})
                    renamed = step_renames.get(out_name, out_name)
                    expr = f'{_call_alias(step_id)}.{renamed}'
                elif binding.get('source') == 'input':
                    bname = binding['name']
                    if input_name:
                        expr = f'{input_name}[{s_var}].{bname}'
                    elif bname in table_columns:
                        expr = f'{_column_input_name(source, bname)}[{s_var}]'
                    else:
                        expr = bname
                elif binding.get('source') == 'literal':
                    expr = json.dumps(binding['value'])
                else:
                    expr = in_name
                lines.append(f'                {in_name} = {expr},')
            lines.append('        }')
    else:
        tname = tasks.get(step.id, _task_name(step.id))
        alias = _call_alias(step.id)
        call_line = f'        call {tname}'
        if alias != tname:
            call_line += f' as {alias}'
        lines.append(call_line + ' {')
        lines.append('            input:')
        task_inputs = (step.task or {}).get('inputs', {})
        for in_name in task_inputs.keys():
            if input_name:
                expr = f'{input_name}[{s_var}].{in_name}'
            else:
                binding = step.bindings.get(in_name)
                if binding is None:
                    if in_name in table_columns:
                        expr = f'{_column_input_name(source, in_name)}[{s_var}]'
                    else:
                        expr = in_name
                elif isinstance(binding, Input):
                    expr = f'{binding.name}[{s_var}]' if binding.name in table_columns else binding.name
                elif isinstance(binding, Literal):
                    expr = _literal_to_wdl(binding.value)
                elif isinstance(binding, Field):
                    base = _binding_to_wdl_expr(binding.source, step.id, None, output_renames)
                    if isinstance(binding.source, StepCall):
                        step_renames = output_renames.get(binding.source.id, {})
                        out_name = step_renames.get(binding.name, binding.name)
                        expr = f'{base}.{out_name}[{s_var}]'
                    else:
                        expr = f'{base}.{binding.name}'
                else:
                    expr = _binding_to_wdl_expr(binding, step.id, None, output_renames)
            lines.append(f'                {in_name} = {expr},')
        lines.append('        }')
    lines.append('    }')
    lines.append('')

    return lines


def _get_table_columns(source):
    return set(table_columns(source).keys())


def _column_input_name(source, col_name):
    return column_input_name(source, col_name)


def _derive_length_expr(source, bindings, table_column_names):
    input_name = source_input_name(source)
    if input_name is not None:
        return f'length({input_name})'
    if source_kind(source) == 'table':
        columns = table_columns(source)
        col_names = list(columns.keys())
        if col_names:
            col = columns[col_names[0]]
            if isinstance(col, Input):
                return f'length({col.name})'
            if isinstance(col, Field) and isinstance(col.source, StepCall):
                return f'length({_call_alias(col.source.id)}.{col.name})'
    for _, binding in bindings.items():
        if isinstance(binding, Input):
            return f'length({binding.name})'
    return '1'


def _mapped_by_step_to_wdl(step, tasks, output_renames=None, inline_map=None):
    if output_renames is None:
        output_renames = {}
    if inline_map is None:
        inline_map = {}
    lines = []
    map_info = step.map or {}
    group_key = map_info.get('group_by')
    source = map_info.get('source', {})
    input_schema = step.input_schema or {}

    col_names = list(input_schema.keys())
    if group_key:
        col_names = [group_key] + [n for n in col_names if n != group_key]

    if step.id in inline_map:
        sub_steps = inline_map[step.id]

        s_var = f'{step.id}_i'
        first_col = next((c for c in col_names), None)
        len_expr = f'length({first_col})' if first_col else '1'
        lines.append(f'    scatter ({s_var} in range({len_expr})) {{')
        for sd in sub_steps:
            sid = sd['id']
            tname = tasks.get(sid, _task_name(sid))
            alias = _call_alias(sid)
            call_line = f'        call {tname}'
            if alias != tname:
                call_line += f' as {alias}'
            lines.append(call_line + ' {')
            lines.append('            input:')
            sub_inputs = sd.get('inputs', {})
            sub_bindings = sd.get('bindings', {})
            for in_name in sub_inputs:
                binding = sub_bindings.get(in_name)
                if binding is None:
                    expr = f'{in_name}[{s_var}]'
                elif binding.get('source') == 'step_output':
                    step_id_ref = binding['step']
                    out_name = binding['output']
                    step_renames = output_renames.get(step_id_ref, {})
                    renamed = step_renames.get(out_name, out_name)
                    expr = f'{_call_alias(step_id_ref)}.{renamed}'
                elif binding.get('source') == 'input':
                    expr = f'{binding["name"]}[{s_var}]'
                elif binding.get('source') == 'literal':
                    expr = json.dumps(binding['value'])
                else:
                    expr = f'{in_name}[{s_var}]'
                lines.append(f'            {in_name} = {expr},')
            lines.append('        }')
        lines.append('    }')
        lines.append('')
        return lines

    grouped_var = f'{step.id}_grouped'

    if source_kind(source) == 'input':
        non_key_names = [n for n in col_names if n != group_key]
        key_t = to_wdl_type(input_schema.get(group_key, 'str'))
        val_t = _val_type(non_key_names, input_schema) if non_key_names else 'String'
        zip_expr = _build_zip_chain_from_names(col_names)
        lines.append(f'    Array[Pair[{key_t}, Array[{val_t}]]] {grouped_var} = collect_by_key({zip_expr})')
    else:
        zip_expr = _build_zip_chain(col_names, source)
        key_t = to_wdl_type(input_schema.get(group_key, 'str'))
        val_t = _val_type(col_names, input_schema)
        lines.append(f'    Array[Pair[{key_t}, Array[{val_t}]]] {grouped_var} = collect_by_key({zip_expr})')

    lines.append('')

    g_var = f'{step.id}_g'
    lines.append(f'    scatter ({g_var} in {grouped_var}) {{')

    lines.append(f'        {to_wdl_type(input_schema.get(group_key, "str"))} {group_key}_val = {g_var}.left')
    for col in col_names:
        if col == group_key:
            continue
        t = to_wdl_type(input_schema.get(col, 'str'))
        path = _col_access_path(col, col_names)
        if path.startswith('right.'):
            path = path[6:]
        lines.append(f'        Array[{t}] {col}_vals = {g_var}.right.{path}')

    tname = tasks.get(step.id, _task_name(step.id))
    alias = _call_alias(step.id)
    call_line = f'        call {tname}'
    if alias != tname:
        call_line += f' as {alias}'
    lines.append(call_line + ' {')
    lines.append('            input:')
    task_inputs = (step.task or {}).get('inputs', {})
    for in_name in task_inputs.keys():
        if in_name == group_key:
            lines.append(f'                {in_name} = {group_key}_val,')
        elif in_name in col_names:
            lines.append(f'                {in_name} = {in_name}_vals,')
        else:
            binding = step.bindings.get(in_name)
            if binding:
                expr = _binding_to_wdl_expr(binding, step.id, None, output_renames)
                lines.append(f'                {in_name} = {expr},')
            else:
                lines.append(f'                {in_name} = {in_name},')
    lines.append('        }')
    lines.append('    }')
    lines.append('')

    return lines


def _build_zip_chain(col_names, source):
    if not col_names:
        return '[]'
    inner = _column_expr(col_names[-1], source)
    for col in reversed(col_names[:-1]):
        col_expr = _column_expr(col, source)
        inner = f'zip({col_expr}, {inner})'
    return inner


def _build_zip_chain_from_names(col_names):
    if not col_names:
        return '[]'
    inner = col_names[-1]
    for col in reversed(col_names[:-1]):
        inner = f'zip({col}, {inner})'
    return inner


def _column_expr(col, source):
    return column_input_name(source, col)


def _val_type(col_names, schema):
    if len(col_names) == 1:
        return to_wdl_type(schema.get(col_names[0], 'str'))
    if len(col_names) == 2:
        return f'Pair[{to_wdl_type(schema.get(col_names[0], "str"))}, {to_wdl_type(schema.get(col_names[1], "str"))}]'
    inner = to_wdl_type(schema.get(col_names[-1], 'str'))
    for col in reversed(col_names[1:-1]):
        inner = f'Pair[{to_wdl_type(schema.get(col, "str"))}, {inner}]'
    return f'Pair[{to_wdl_type(schema.get(col_names[0], "str"))}, {inner}]'


def _col_access_path(col, col_names):
    idx = col_names.index(col)
    if idx == 0:
        return 'left'
    path = 'right'
    for _ in range(1, idx):
        path += '.right'
    path += '.left'
    return path


def _subworkflow_to_wdl(step, parent_id):
    definition = step.task or {}
    dag_data = definition.get('dag', {})
    if not dag_data:
        return ''
    wf_id = f'{parent_id}_{step.id}'
    return transpile_dag_dict(dag_data, workflow_id=wf_id, _top_level=False)


def _collect_structs(dag):
    structs = []
    seen = set()
    for step in dag.steps:
        for binding in step.bindings.values():
            if isinstance(binding, Record):
                shape = tuple(sorted(binding.fields.keys()))
                if shape and shape not in seen:
                    seen.add(shape)
                    structs.append(_emit_struct(binding))
    return structs


def _struct_name_for(binding):
    keys = tuple(sorted(binding.fields.keys()))
    return '_Rec_' + '_'.join(k.capitalize() for k in keys)


def _infer_fieldto_wdl_type(fbinding):
    if isinstance(fbinding, Input):
        return to_wdl_type(fbinding.type or 'str')
    if isinstance(fbinding, Literal):
        return _infer_literal_type(fbinding.value)
    if isinstance(fbinding, Field):
        if isinstance(fbinding.source, StepCall):
            out_type = (fbinding.source.task or {}).get('outputs', {}).get(fbinding.name, {}).get('type', 'str')
            return to_wdl_type(out_type)
        if isinstance(fbinding.source, Input):
            return to_wdl_type(fbinding.source.type or 'str')
    return 'String'


def _emit_struct(record):
    name = _struct_name_for(record)
    lines = [f'struct {name} {{']
    for fname in sorted(record.fields.keys()):
        typ = _infer_fieldto_wdl_type(record.fields[fname])
        lines.append(f'    {typ} {fname}')
    lines.append('}')
    return '\n'.join(lines)


def _input_struct_name(step_id):
    name = step_name(step_id, default='Input')
    return name[0].upper() + name[1:] + 'Input'


def _emit_input_struct(step_id, input_schema):
    name = _input_struct_name(step_id)
    lines = [f'struct {name} {{']
    for col_name in sorted(input_schema):
        wdl_type = to_wdl_type(input_schema[col_name])
        lines.append(f'    {wdl_type} {col_name}')
    lines.append('}')
    return '\n'.join(lines)


def _infer_literal_type(value):
    if isinstance(value, bool):
        return 'Boolean'
    if isinstance(value, int):
        return 'Int'
    if isinstance(value, float):
        return 'Float'
    if value is None:
        return 'String?'
    return 'String'
