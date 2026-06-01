import json

from swl.ir.dag import DAG, Field, Input, Literal, Merge, Record, StepCall


def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())


def transpile_dag_dict(data, workflow_id='main', _top_level=True):
    dag = DAG.from_dict(data)
    _validate_supported(dag)

    structs = _collect_structs(dag)
    tasks = {}
    sub_workflows = []
    for step in dag.steps:
        if step.id not in tasks:
            if step.type == 'workflow':
                wdl = _subworkflow_to_wdl(step, workflow_id)
                tasks[step.id] = _task_name(step.id)
                sub_workflows.append(wdl)
            else:
                tasks[step.id] = _task_name(step.id)

    lines = []
    if _top_level:
        lines.append('version 1.1\n')
    for s in structs:
        lines.append(s)
        lines.append('')
    for step in dag.steps:
        if step.type != 'workflow' and step.id not in [s.id for s in dag.steps if s.type == 'workflow']:
            lines.append(_task_to_wdl(step))
            lines.append('')
    for sw in sub_workflows:
        lines.append(sw)
        lines.append('')
    lines.append(_dag_to_wdl(dag, workflow_id, tasks))
    return '\n'.join(lines)


def _wdl_type(swl_type, optional=False):
    base = {
        'file': 'File',
        'str': 'String',
        'int': 'Int',
        'float': 'Float',
        '[file]': 'Array[File]',
        '[str]': 'Array[String]',
        '[int]': 'Array[Int]',
        '[float]': 'Array[Float]',
    }.get(swl_type, 'String')
    if optional:
        return base + '?'
    return base


def _task_name(step_id):
    name = step_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'task'
    if name[0].isdigit():
        name = '_' + name
    return name.lower()


def _wf_name(workflow_id):
    name = workflow_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'main'
    if name[0].isdigit():
        name = '_' + name
    return name


def _call_alias(step_id):
    name = step_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'step'
    if name[0].isdigit():
        name = '_' + name
    return name


def _task_to_wdl(step):
    task = step.task or {}
    body = task.get('body', '')
    lines = [f'task {_task_name(step.id)} {{', '']

    inputs = task.get('inputs', {})
    if inputs:
        lines.append('    input {')
        for name, spec in inputs.items():
            t = _wdl_type(spec.get('type'))
            lines.append(f'        {t} {name}')
        lines.append('    }')
        lines.append('')

    lines.append('    command <<<')
    for line in body.split('\n'):
        interp_line = _interpolate_bash_vars(line)
        lines.append(f'        {interp_line}')
    lines.append('    >>>')
    lines.append('')

    outputs = task.get('outputs', {})
    if outputs:
        lines.append('    output {')
        for name, spec in outputs.items():
            t = _wdl_type(spec.get('type'))
            default = spec.get('default')
            if default:
                path_expr = _interp_to_wdl(default)
                lines.append(f'        {t} {name} = "{path_expr}"')
            else:
                lines.append(f'        {t} {name} = "*.{name}"')
        lines.append('    }')
        lines.append('')

    req = _emit_requirements(step)
    if req:
        lines.append(req)
        lines.append('')

    lines.append('}')
    return '\n'.join(lines)


def _interpolate_bash_vars(line):
    import re
    result = []
    i = 0
    while i < len(line):
        if line[i] == '$' and i + 1 < len(line) and line[i + 1] == '{':
            brace_start = i
            depth = 0
            for j in range(i, len(line)):
                if line[j] == '{':
                    depth += 1
                elif line[j] == '}':
                    depth -= 1
                    if depth == 0:
                        var_name = line[brace_start + 2:j]
                        result.append(f'~{{{var_name}}}')
                        i = j + 1
                        break
            else:
                result.append(line[i])
                i += 1
        else:
            result.append(line[i])
            i += 1
    return ''.join(result)


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
    if value is None:
        return None
    if value.get('kind') == 'word':
        parts = value.get('parts', [])
        result = ''
        for part in parts:
            if part.get('kind') == 'literal':
                result += part['text']
            elif part.get('kind') == 'var':
                result += f"~{{{part['name']}}}"
            elif part.get('kind') == 'expr':
                result += f"~{{{part['text']}}}"
        return result
    return None


def _dag_to_wdl(dag, workflow_id, tasks):
    lines = [f'workflow {_wf_name(workflow_id)} {{', '']

    if dag.inputs:
        lines.append('    input {')
        for name, spec in dag.inputs.items():
            t = _wdl_type(spec.type)
            lines.append(f'        {t} {name}')
        lines.append('    }')
        lines.append('')

    for step in dag.steps:
        if getattr(step, 'map', None) is not None:
            if step.map.get('group_by') is not None:
                lines.extend(_mapped_by_step_to_wdl(step, tasks))
            else:
                lines.extend(_mapped_step_to_wdl(step, tasks))
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
            expr = _binding_to_wdl_expr(binding, step.id, dag)
            lines.append(f'            {in_name} = {expr},')
        lines.append('    }')
        lines.append('')

    if dag.outputs:
        lines.append('    output {')
        for name, binding in dag.outputs.items():
            expr = _binding_to_wdl_expr(binding, None, dag)
            typ = _infer_output_type(name, binding, dag)
            lines.append(f'        {typ} {name} = {expr}')
        lines.append('    }')

    lines.append('}')
    return '\n'.join(lines)


def _binding_to_wdl_expr(binding, current_step_id, dag):
    if isinstance(binding, Input):
        return binding.name

    if isinstance(binding, Literal):
        return _literal_to_wdl(binding.value)

    if isinstance(binding, Field):
        source_expr = _binding_to_wdl_expr(binding.source, current_step_id, dag)
        return f'{source_expr}.{binding.name}'

    if isinstance(binding, Merge):
        raise ValueError(
            f'Merge bindings must be flattened before WDL transpilation'
        )

    if isinstance(binding, Record):
        struct_name = _struct_name_for(binding)
        fields = []
        for fname, fbinding in binding.fields.items():
            fexpr = _binding_to_wdl_expr(fbinding, current_step_id, dag)
            fields.append(f'{fname}: {fexpr}')
        return f'{struct_name} {{{", ".join(fields)}}}'

    if isinstance(binding, StepCall):
        return _call_alias(binding.id)

    if isinstance(binding, dict):
        return _dict_binding_to_wdl_expr(binding, current_step_id, dag)

    raise ValueError(f'Unsupported binding for WDL: {type(binding).__name__}')


def _dict_binding_to_wdl_expr(binding, current_step_id, dag):
    source = binding.get('source')
    if source == 'input':
        return binding['name']
    if source == 'literal':
        return _literal_to_wdl(binding.get('value'))
    if source == 'field':
        inner = _dict_binding_to_wdl_expr(binding['value'], current_step_id, dag)
        return f'{inner}.{binding["field"]}'
    if source == 'step_call':
        return _call_alias(binding['step'])
    if 'step' in binding and 'output' in binding:
        alias = _call_alias(binding['step'])
        return f'{alias}.{binding["output"]}'
    if source == 'record':
        struct_name = _record_struct_name(binding)
        fields = []
        for fname, fbinding in binding.get('fields', {}).items():
            fexpr = _dict_binding_to_wdl_expr(fbinding, current_step_id, dag)
            fields.append(f'{fname}: {fexpr}')
        return f'{struct_name} {{{", ".join(fields)}}}'
    if source == 'merge':
        raise ValueError(
            f'Merge bindings must be flattened before WDL transpilation'
        )
    raise ValueError(f'Unsupported binding source for WDL: {source!r}')


def _literal_to_wdl(value):
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if value is None:
        return 'None'
    return str(value)


def _record_struct_name(binding):
    fields = binding.get('fields', {})
    keys = tuple(sorted(fields.keys()))
    return '_Rec_' + '_'.join(k.capitalize() for k in keys)


def _mapped_step_to_wdl(step, tasks):
    lines = []
    tname = tasks.get(step.id, _task_name(step.id))
    alias = _call_alias(step.id)
    s_var = f'{step.id}_i'

    map_info = step.map or {}
    source = map_info.get('source', {})
    table_columns = _get_table_columns(source)
    len_expr = _derive_length_expr(source, step.bindings, table_columns)
    lines.append(f'    scatter ({s_var} in range({len_expr})) {{')

    call_line = f'        call {tname}'
    if alias != tname:
        call_line += f' as {alias}'
    lines.append(call_line + ' {')
    lines.append('            input:')
    task_inputs = (step.task or {}).get('inputs', {})
    for in_name in task_inputs.keys():
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
            base = _binding_to_wdl_expr(binding.source, step.id, None)
            if isinstance(binding.source, StepCall):
                expr = f'{base}.{binding.name}[{s_var}]'
            else:
                expr = f'{base}.{binding.name}'
        else:
            expr = _binding_to_wdl_expr(binding, step.id, None)
        lines.append(f'                {in_name} = {expr},')
    lines.append('        }')
    lines.append('    }')
    lines.append('')

    return lines


def _get_table_columns(source):
    if isinstance(source, dict):
        if source.get('source') == 'table':
            return set(source.get('columns', {}).keys())
    return set()


def _column_input_name(source, col_name):
    if isinstance(source, dict):
        if source.get('source') == 'input':
            return source['name']
        if source.get('source') == 'table':
            col = source.get('columns', {}).get(col_name, {})
            if isinstance(col, dict) and col.get('source') == 'input':
                return col['name']
    return col_name


def _derive_length_expr(source, bindings, table_columns):
    if isinstance(source, dict):
        if source.get('source') == 'input':
            return f'length({source["name"]})'
        if source.get('source') == 'table':
            columns = source.get('columns', {})
            col_names = list(columns.keys())
            if col_names:
                col = columns[col_names[0]]
                if isinstance(col, Input):
                    return f'length({col.name})'
                if isinstance(col, Field) and isinstance(col.source, StepCall):
                    return f'length({_call_alias(col.source.id)}.{col.name})'
    for name, binding in bindings.items():
        if isinstance(binding, Input):
            return f'length({binding.name})'
    return '1'


def _mapped_by_step_to_wdl(step, tasks):
    lines = []
    tname = tasks.get(step.id, _task_name(step.id))
    alias = _call_alias(step.id)
    map_info = step.map or {}
    group_key = map_info.get('group_by')
    source = map_info.get('source', {})
    input_schema = step.input_schema or {}

    col_names = list(input_schema.keys())
    if group_key and group_key not in col_names:
        col_names = [group_key] + col_names

    zip_expr = _build_zip_chain(col_names, source)
    grouped_var = f'{step.id}_grouped'

    key_t = _wdl_type(input_schema.get(group_key, 'str'))
    val_t = _val_type(col_names, input_schema)
    lines.append(f'    Array[Pair[{key_t}, Array[{val_t}]]] {grouped_var} = collect_by_key({zip_expr})')
    lines.append('')

    g_var = f'{step.id}_g'
    lines.append(f'    scatter ({g_var} in {grouped_var}) {{')

    lines.append(f'        {_wdl_type(input_schema.get(group_key, "str"))} {group_key}_val = {g_var}.left')
    for col in col_names:
        if col == group_key:
            continue
        t = _wdl_type(input_schema.get(col, 'str'))
        lines.append(f'        Array[{t}] {col}_vals = {g_var}.right.{_col_access_path(col, col_names)}')

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
                expr = _binding_to_wdl_expr(binding, step.id, None)
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


def _column_expr(col, source):
    if isinstance(source, dict):
        if source.get('source') == 'input':
            return source['name']
        columns = source.get('columns', {})
        col_binding = columns.get(col, {})
        if isinstance(col_binding, dict) and col_binding.get('source') == 'input':
            return col_binding['name']
    return col


def _key_type(key, schema):
    return _wdl_type(schema.get(key, 'str'))


def _val_type(col_names, schema):
    if len(col_names) == 1:
        return _wdl_type(schema.get(col_names[0], 'str'))
    if len(col_names) == 2:
        return f'Pair[{_wdl_type(schema.get(col_names[0], "str"))}, {_wdl_type(schema.get(col_names[1], "str"))}]'
    inner = _wdl_type(schema.get(col_names[-1], 'str'))
    for col in reversed(col_names[1:-1]):
        inner = f'Pair[{_wdl_type(schema.get(col, "str"))}, {inner}]'
    return f'Pair[{_wdl_type(schema.get(col_names[0], "str"))}, {inner}]'


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
                shape = _record_shape(binding)
                if shape and shape not in seen:
                    seen.add(shape)
                    structs.append(_emit_struct(shape))
            if isinstance(binding, dict) and binding.get('source') == 'record':
                shape = _dict_record_shape(binding)
                if shape and shape not in seen:
                    seen.add(shape)
                    structs.append(_emit_struct(shape))
    return structs


def _record_shape(binding):
    if not isinstance(binding, Record):
        return None
    keys = tuple(sorted(binding.fields.keys()))
    return keys


def _dict_record_shape(binding):
    if not isinstance(binding, dict) or binding.get('source') != 'record':
        return None
    keys = tuple(sorted(binding.get('fields', {}).keys()))
    return keys


def _struct_name_for(binding):
    if isinstance(binding, Record):
        keys = tuple(sorted(binding.fields.keys()))
    elif isinstance(binding, dict):
        keys = tuple(sorted(binding.get('fields', {}).keys()))
    else:
        keys = binding
    return '_Rec_' + '_'.join(k.capitalize() for k in keys)


def _emit_struct(shape):
    name = _struct_name_for(shape)
    lines = [f'struct {name} {{']
    for field_name in shape:
        lines.append(f'    String {field_name}')
    lines.append('}')
    return '\n'.join(lines)


def _infer_output_type(name, binding, dag):
    if isinstance(binding, Field) and isinstance(binding.source, StepCall):
        step = binding.source
        out_type = (step.task or {}).get('outputs', {}).get(binding.name, {}).get('type', 'str')
        wdl_t = _wdl_type(out_type)
        if getattr(step, 'map', None) is not None:
            return f'Array[{wdl_t}]'
        return wdl_t
    if isinstance(binding, Field) and isinstance(binding.source, Input):
        input_spec = dag.inputs.get(binding.source.name)
        if input_spec:
            return _wdl_type(input_spec.type)
        return 'String'
    if isinstance(binding, Input):
        spec = dag.inputs.get(binding.name)
        return _wdl_type(spec.type if spec else 'str')
    if isinstance(binding, Literal):
        return _infer_literal_type(binding.value)
    if isinstance(binding, dict):
        return _dict_infer_output_type(binding, dag)
    return 'String'


def _dict_infer_output_type(binding, dag):
    source = binding.get('source')
    if source == 'input':
        spec = dag.inputs.get(binding['name'])
        return _wdl_type(spec.type if spec else 'str')
    if source == 'literal':
        return _infer_literal_type(binding.get('value'))
    if 'step' in binding and 'output' in binding:
        return 'String'
    return 'String'


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


def _validate_supported(dag):
    for step in dag.steps:
        for name, value in step.bindings.items():
            _validate_binding(value, step.id)

    for name, value in dag.outputs.items():
        _validate_output_binding(value, name)


def _validate_binding(value, step_id):
    if isinstance(value, Merge):
        raise ValueError(
            f'WDL transpilation does not support merge bindings '
            f'(step {step_id}). Flatten merges before transpiling.'
        )
    if isinstance(value, dict) and value.get('source') == 'merge':
        raise ValueError(
            f'WDL transpilation does not support merge bindings '
            f'(step {step_id}). Flatten merges before transpiling.'
        )


def _validate_output_binding(value, name):
    if isinstance(value, Literal):
        raise ValueError(f'WDL does not support literal workflow outputs: {name}')
    if isinstance(value, Merge):
        raise ValueError(f'WDL does not support merge workflow outputs: {name}')
