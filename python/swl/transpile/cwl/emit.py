import json
import re

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import (
    classify_var,
    field_chain_parts,
    field_path_after_first,
    flatten_dag_outputs,
    run_value,
    source_input_name,
    source_kind,
    table_columns,
    validate_dag_for_transpile,
    word_interp,
)
from swl.types import to_array_type, to_cwl_type


def transpile_dag_file(path):
    data = json.load(open(path))
    return transpile_dag_dict(data, workflow_id='main')


def transpile_dag_dict(data, workflow_id='main'):
    dag = DAG.from_dict(data)
    dag.validate()
    flatten_dag_outputs(dag)
    _validate_supported(dag)
    tools = []
    tool_ids = {}
    map_by_tools = {}
    map_by_steps = {}

    record_tools = []
    record_map = {}
    record_output_map = {}

    for step in dag.steps:
        for name, value in list(step.bindings.items()):
            if isinstance(value, Record):
                tool, step_entry, source = _emit_record_tool(workflow_id, step.id, name, value, dag)
                record_tools.append(tool)
                record_map[(step.id, name)] = (step_entry, source)

    for name, output in list(dag.outputs.items()):
        value = output.value if isinstance(output, OutputSpec) else output
        if isinstance(value, Record):
            tool, step_entry, source = _emit_record_tool(workflow_id, 'outputs', name, value, dag)
            record_tools.append(tool)
            record_output_map[name] = (step_entry, source)

    for step in dag.steps:
        tool_id = step.id
        if getattr(step, 'map', None) is not None and step.map.get('group_by') is not None:
            map_graph, map_steps = _emit_map_by_graph(step, dag)
            map_by_tools[tool_id] = map_graph
            map_by_steps[tool_id] = map_steps
            continue
        if tool_id not in tool_ids:
            tool_ids[tool_id] = f'#{tool_id}'
            tools.extend(_tool_to_cwl(step, tool_ids[tool_id]))

    workflow_inputs = dict(dag.inputs)
    for step in dag.steps:
        if getattr(step, 'map', None) is None:
            continue
        source = step.map.get('source', {})
        schema = (getattr(step, 'input_schema', None) or {})
        for name, typ in schema.items():
            if name not in workflow_inputs:
                workflow_inputs[name] = {'type': to_array_type(typ), 'desc': None}
            elif source_kind(source) == 'table':
                column = table_columns(source).get(name, {})
                if isinstance(column, dict) and column.get('source') == 'input' and column.get('name') == name:
                    workflow_inputs[name] = {'type': to_array_type(typ), 'desc': None}
        if source_input_name(source) is not None:
            src_name = source_input_name(source)
            if src_name in workflow_inputs:
                for name, typ in schema.items():
                    workflow_inputs[name] = {'type': to_array_type(typ), 'desc': None}

    all_tools = tools[:]
    for mt in map_by_tools.values():
        all_tools.extend(mt)
    all_tools.extend(record_tools)

    def step_to_cwl(step):
        if step.id in map_by_steps:
            return map_by_steps[step.id]
        steps_out = []
        for (s_id, b_name), (step_entry, _) in record_map.items():
            if s_id == step.id:
                steps_out.append(step_entry)
        steps_out.append(_step_to_cwl(workflow_id, step, tool_ids[step.id], record_map))
        return steps_out

    wf_steps = []
    for step in dag.steps:
        wf_steps.extend(step_to_cwl(step))
    for _, (step_entry, _) in record_output_map.items():
        wf_steps.insert(0, step_entry)

    def _is_map_by_step_ref(value, dag):
        if isinstance(value, Field) and isinstance(value.source, StepCall) and getattr(value.source, 'map', None):
            return bool(value.source.map.get('group_by'))
        return False

    outputs = []
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if name in record_output_map:
            _, source = record_output_map[name]
            outputs.append({
                'id': f'#{workflow_id}/{name}',
                'type': to_cwl_type(output.type),
                'outputSource': source,
            })
        else:
            source = _binding_source(value, workflow_id)
            typ = to_cwl_type(output.type)
            if _is_map_by_step_ref(value, dag):
                typ = {'type': 'array', 'items': typ}
            outputs.append({
                'id': f'#{workflow_id}/{name}',
                'type': typ,
                'outputSource': source,
            })

    workflow = {
        'id': f'#{workflow_id}',
        'class': 'Workflow',
        'inputs': [_workflow_input_to_cwl(workflow_id, name, spec) for name, spec in workflow_inputs.items()],
        'outputs': outputs,
        'requirements': [{'class': 'ScatterFeatureRequirement'}, {'class': 'MultipleInputFeatureRequirement'}, {'class': 'SubworkflowFeatureRequirement'}, {'class': 'StepInputExpressionRequirement'}],
        'steps': wf_steps,
    }
    return {
        'cwlVersion': 'v1.2',
        '$graph': all_tools + [workflow],
    }


def _tool_to_cwl(step, tool_id):
    definition = step.task or {}
    if definition.get('class') == 'Workflow':
        packed = transpile_dag_dict(definition['dag'], workflow_id=tool_id[1:])
        return packed.get('$graph', [])
    run = definition.get('run', {})
    body = definition.get('body', '')
    entry = _interpolate_shell(body, step)
    requirements = [
        {
            'class': 'InitialWorkDirRequirement',
            'listing': [
                {
                    'entryname': 'script.sh',
                    'entry': entry,
                }
            ],
        }
    ]
    resource = _resource_requirement(run)
    if resource is not None:
        requirements.append(resource)
    docker = _docker_requirement(run)
    if docker is not None:
        requirements.append(docker)
    hints = _hints_from_run(run)
    needs_js = '$' in body
    if not needs_js:
        needs_js = any(_has_expr_interpolation(spec) for spec in definition.get('outputs', {}).values())
    if needs_js:
        requirements.append({'class': 'InlineJavascriptRequirement'})
    input_names = list(definition.get('inputs', {}).keys())
    tool = {
        'id': tool_id,
        'class': 'CommandLineTool',
        'baseCommand': ['bash', 'script.sh'],
        'inputs': [_tool_input_to_cwl(tool_id, name, spec) for name, spec in definition.get('inputs', {}).items()],
        'outputs': [_tool_output_to_cwl(tool_id, name, spec, input_names) for name, spec in definition.get('outputs', {}).items()],
        'requirements': requirements,
    }
    if hints:
        tool['hints'] = hints
    return [tool]


def _workflow_input_to_cwl(workflow_id, name, spec):
    typ = spec.type if hasattr(spec, 'type') else spec.get('type')
    desc = spec.desc if hasattr(spec, 'desc') else spec.get('desc')
    return {
        'id': f'#{workflow_id}/{name}',
        'type': to_cwl_type(typ),
        **({'doc': desc} if desc else {}),
    }


def _step_to_cwl(workflow_id, step, tool_id, record_map=None):
    if record_map is None:
        record_map = {}
    inputs = []
    for name, binding in step.bindings.items():
        if (step.id, name) in record_map:
            _, source = record_map[(step.id, name)]
            inputs.append({
                'id': f'#{workflow_id}/{step.id}/{name}',
                'source': source,
            })
        else:
            inputs.append(_step_input_to_cwl(workflow_id, step.id, name, binding))
    data = {
        'id': f'#{workflow_id}/{step.id}',
        'run': tool_id,
        'in': inputs,
        'out': [f'#{workflow_id}/{step.id}/{name}' for name in step.outputs],
    }
    if getattr(step, 'map', None) is not None:
        scatter_ports = step.map.get('scatter') or sorted((getattr(step, 'input_schema', None) or {}).keys())
        if scatter_ports:
            source = step.map.get('source', {})
            if source_input_name(source) is not None:
                for port in scatter_ports:
                    if not any(item['id'] == f'#{workflow_id}/{step.id}/{port}' for item in data['in']):
                        data['in'].append({'id': f'#{workflow_id}/{step.id}/{port}', 'source': f'#{workflow_id}/{port}'})
            elif source_kind(source) == 'table':
                columns = table_columns(source)
                for port in scatter_ports:
                    column = columns.get(port)
                    source_ref = f'#{workflow_id}/{port}'
                    if isinstance(column, dict) and source_kind(column) == 'input' and 'name' in column:
                        source_ref = f"#{workflow_id}/{column['name']}"
                    if not any(item['id'] == f'#{workflow_id}/{step.id}/{port}' for item in data['in']):
                        data['in'].append({'id': f'#{workflow_id}/{step.id}/{port}', 'source': source_ref})
            data['scatter'] = [f'#{workflow_id}/{step.id}/{port}' for port in scatter_ports]
            data['scatterMethod'] = 'dotproduct'
    return data


def _step_input_to_cwl(workflow_id, task_id, name, value):
    if isinstance(value, Literal):
        return {
            'id': f'#{workflow_id}/{task_id}/{name}',
            'default': value.value,
        }
    result = {
        'id': f'#{workflow_id}/{task_id}/{name}',
        'source': _binding_source(value, workflow_id),
    }
    if isinstance(value, Field) and value.name != name:
        result['valueFrom'] = f"$({value.name})"
    return result


def _tool_input_to_cwl(tool_id, name, spec):
    return {
        'id': f'{tool_id}/{name}',
        'type': to_cwl_type(spec.get('type')),
        **({'doc': spec.get('desc')} if spec.get('desc') else {}),
    }


def _tool_output_to_cwl(tool_id, name, spec, input_names=()):
    return {
        'id': f'{tool_id}/{name}',
        'type': to_cwl_type(spec.get('type')),
        'outputBinding': {
            'glob': _interp_to_cwl_glob(spec.get('default'), input_names),
        },
        **({'doc': spec.get('desc')} if spec.get('desc') else {}),
    }


def _resource_requirement(run):
    req = {'class': 'ResourceRequirement'}
    cpu = run_value(run, 'cpu')
    memory = run_value(run, 'memory')
    if cpu is not None:
        req['coresMin'] = cpu
    if memory is not None:
        req['ramMin'] = memory
    if len(req) == 1:
        return None
    return req


def _hints_from_run(run):
    time = run_value(run, 'time')
    return [{'class': 'TimeLimit', 'timeLimit': time}] if time is not None else []


def _docker_requirement(run):
    image = run_value(run, 'image')
    if image is None:
        return None
    return {'class': 'DockerRequirement', 'dockerPull': image}


def _binding_source(value, workflow_id='main'):
    if isinstance(value, Input):
        return f'#{workflow_id}/{value.name}'
    if isinstance(value, Literal):
        return value.value
    if isinstance(value, StepCall):
        return f'#{workflow_id}/{value.id}'
    if isinstance(value, Field):
        src = value.source
        if isinstance(src, Input):
            return f'#{workflow_id}/{src.name}/{value.name}'
        if isinstance(src, StepCall):
            return f'#{workflow_id}/{src.id}/{value.name}'
        if isinstance(src, Field):
            root, root_name, _ = field_chain_parts(value)
            if isinstance(root, Input):
                return f'#{workflow_id}/{root.name}'
            if isinstance(root, StepCall):
                return f'#{workflow_id}/{root.id}/{root_name}'
    raise ValueError(f'Unsupported value for CWL binding source: {type(value).__name__} {value!r}')


def _validate_supported(dag):
    for step in dag.steps:
        if getattr(step, 'map', None) is not None and step.map.get('group_by') is not None:
            _validate_map_by_preconditions(step)
    validate_dag_for_transpile(dag, 'CWL')
    for step in dag.steps:
        for name, spec in step.task.get('outputs', {}).items():
            try:
                _interp_to_cwl_glob(spec.get('default'), step.task.get('inputs', {}).keys())
            except ValueError as exc:
                raise ValueError(f'Unsupported step output path for CWL transpilation: {step.id}.{name}: {exc}') from exc
    for step in dag.steps:
        for binding in step.bindings.values():
            if isinstance(binding, Merge):
                raise ValueError(
                    f'CWL does not support Merge bindings: step {step.id} contains a Merge value'
                )


def _interp_to_cwl_glob(value, input_names=()):
    if value is None:
        return '*'
    input_set = set(input_names) if input_names else set()

    def _resolve_expr(text):
        if not input_set:
            return f"({text})"
        resolved = text
        for var in re.findall(r'[A-Za-z_]\w*', text):
            if var in input_set:
                resolved = resolved.replace(var, f'inputs.{var}', 1)
        return f"({resolved})"

    rendered = word_interp(value, lambda text: repr(text), lambda name: f"inputs.{name}", _resolve_expr, joiner=' + ')
    if rendered is None:
        raise ValueError(f'Unsupported interpolation for CWL glob: {value!r}')
    return '$(' + rendered + ')'


def _has_expr_interpolation(spec):
    default = spec.get('default') if isinstance(spec, dict) else getattr(spec, 'default', None)
    if default is None:
        return False
    if isinstance(default, dict) and default.get('kind') == 'word':
        return any(part.get('kind') == 'expr' for part in default.get('parts', []))
    return False


def _interpolate_shell(body, step):
    if '$' not in body:
        return body
    task = step.task or {}
    input_names = set(task.get('inputs', {}).keys())
    input_types = {n: s.get('type') for n, s in task.get('inputs', {}).items()}
    run = task.get('run', {})
    run_values = {}
    for rv in ('cpu', 'memory', 'time'):
        spec = run.get(rv, {})
        if isinstance(spec, dict) and spec.get('value') is not None:
            run_values[rv] = spec['value']
    run_names = set(run_values.keys())
    output_names = set(step.outputs)

    def _resolve(name):
        scope = classify_var(name, input_names, output_names, run_names)
        if scope == 'input':
            typ = input_types.get(name)
            if typ == '[file]':
                return f'inputs.{name}.map(function(f){{return f.path;}}).join(" ")'
            if typ == 'file':
                return f'inputs.{name}.path'
            return f'inputs.{name}'
        if scope == 'run':
            val = run_values.get(name)
            if val is not None:
                return json.dumps(val)
        return None

    def _line_to_js(line):
        if '$' not in line:
            return json.dumps(line)
        parts = []
        last_end = 0
        for m in re.finditer(r'(?<!\$)\$\{(.+?)\}|(?<!\$)\$(\w+)', line):
            if m.start() > last_end:
                parts.append(json.dumps(line[last_end:m.start()]))
            if m.group(1) is not None:
                content = m.group(1).strip()
                if re.fullmatch(r'\w+', content):
                    r = _resolve(content)
                    parts.append(json.dumps(m.group(0)) if r is None else r)
                else:
                    resolved = content
                    for var in re.findall(r'[A-Za-z_]\w*', content):
                        r = _resolve(var)
                        if r is not None:
                            resolved = resolved.replace(var, r, 1)
                    parts.append(f'({resolved})')
            else:
                r = _resolve(m.group(2))
                parts.append(json.dumps(m.group(0)) if r is None else r)
            last_end = m.end()
        if last_end < len(line):
            parts.append(json.dumps(line[last_end:]))
        if not parts:
            return json.dumps(line)
        return ' + '.join(parts)

    js_lines = [_line_to_js(line) for line in body.split('\n')]
    if len(js_lines) == 1:
        return '$(' + js_lines[0] + ')'
    return '$(' + ' + "\\n" + '.join(js_lines) + ')'


# record bindings ---------------------------------------------------------

def _infer_record_field_type(field_value, dag):
    if isinstance(field_value, Input):
        spec = dag.inputs.get(field_value.name)
        typ = spec.type if hasattr(spec, 'type') else None
        return to_cwl_type(typ)
    if isinstance(field_value, Literal):
        return to_cwl_type(type(field_value.value).__name__)
    if isinstance(field_value, (Field, StepCall)):
        return 'string'
    return 'string'


def _emit_record_tool(workflow_id, step_id, binding_name, record, dag):
    tool_id = f'#rec_{step_id}_{binding_name}'
    step_entry_id = f'#{workflow_id}/rec_{step_id}_{binding_name}'
    field_names = sorted(record.fields.keys())
    tool_inputs = []
    step_inputs = []
    for fname in field_names:
        fvalue = record.fields[fname]
        typ = _infer_record_field_type(fvalue, dag)
        tool_inputs.append({'id': f'{tool_id}/{fname}', 'type': typ})
        if isinstance(fvalue, Literal):
            step_inputs.append({'id': f'{step_entry_id}/{fname}', 'default': fvalue.value})
        else:
            step_inputs.append({'id': f'{step_entry_id}/{fname}', 'source': _binding_source(fvalue, workflow_id)})
    js_lines = ['var fs = require("fs");']
    for fname in field_names:
        js_lines.append(f'record.{fname} = inputs.{fname};')
    js_lines.append('fs.writeFileSync("record.json", JSON.stringify(record));')
    js_lines.append('return {"class": "File", "path": "record.json"};')
    js_expression = '${' + ('var record = {};\n' + '\n'.join(js_lines)) + '}'
    tool = {
        'id': tool_id,
        'class': 'ExpressionTool',
        'requirements': [{'class': 'InlineJavascriptRequirement'}],
        'inputs': tool_inputs,
        'outputs': [{
            'id': f'{tool_id}/record_out',
            'type': 'File',
        }],
        'expression': js_expression,
    }
    step_entry = {
        'id': step_entry_id,
        'run': tool_id,
        'in': step_inputs,
        'out': [f'{step_entry_id}/record_out'],
    }
    return tool, step_entry, f'{step_entry_id}/record_out'


# map_by (grouped scatter) ------------------------------------------------

def _validate_map_by_preconditions(step):
    map_info = step.map or {}
    group_by = map_info.get('group_by')
    schema = step.input_schema or {}
    if group_by not in schema and group_by not in (map_info.get('source', {}).get('columns', {}) or {}):
        raise ValueError(
            f'map_by key {group_by!r} must be a column in the input schema'
        )
    for name, typ in schema.items():
        if typ.startswith('['):
            raise ValueError(
                f'map_by does not support array-typed input {name!r} ({typ}) in mapped function'
            )
    for name, spec in (getattr(step, 'task', None) or {}).get('outputs', {}).items():
        if _has_expr_interpolation(spec):
            raise ValueError(
                f'map_by does not support expr interpolation in mapped function output {name!r}'
            )


def _topological_sort(steps):
    deps = {s['id']: set(s.get('deps', [])) for s in steps}
    ordered = []
    while deps:
        ready = [sid for sid, d in deps.items() if not d]
        if not ready:
            break
        ready.sort()
        for sid in ready:
            ordered.append(sid)
            del deps[sid]
            for d in deps.values():
                d.discard(sid)
    remaining = list(deps.keys())
    remaining.sort()
    ordered.extend(remaining)
    id_map = {s['id']: s for s in steps}
    return [id_map[sid] for sid in ordered]


def _default_to_bash_expr(default):
    if not isinstance(default, dict):
        return None
    if default.get('kind') == 'word':
        parts = default.get('parts', [])
        result = ''
        for part in parts:
            if part.get('kind') == 'var':
                result += '"${' + part['name'] + '}"'
            elif part.get('kind') == 'literal':
                result += part.get('text', '')
            else:
                result += str(part)
        return result
    return None


def _dag_to_combined_script(dag_data):
    steps = list(dag_data.get('steps', []))
    if not steps:
        return ''
    ordered = _topological_sort(steps)
    body_lines = []
    for step in ordered:
        script = step.get('script', '')
        if script:
            if body_lines:
                body_lines.append('')
            body_lines.append(script)
        for oname, ospec in step.get('outputs', {}).items():
            default = ospec.get('default', {})
            bash_expr = _default_to_bash_expr(default)
            if bash_expr:
                body_lines.append(f'export {oname}={bash_expr}')
    return '\n'.join(body_lines)


def _emit_map_by_graph(step, dag):
    map_info = step.map or {}
    group_by = map_info.get('group_by')
    schema = step.input_schema or {}
    output_schema = step.output_schema or {}
    step_id = step.id
    group_tool_id = f'#group_{step_id}'
    wrapper_tool_id = f'#wrap_{step_id}'

    col_names = sorted(schema.keys())
    key_var = 'group_by_key'

    grouping_inputs = [{
        'id': f'{group_tool_id}/{key_var}',
        'type': 'string',
    }]
    for col in col_names:
        grouping_inputs.append({
            'id': f'{group_tool_id}/{col}',
            'type': {'type': 'array', 'items': to_cwl_type(schema[col])},
        })

    listing = [{'entryname': 'inputs.json', 'entry': '$(JSON.stringify(inputs))'}]
    group_script = (
        'import json, os\n'
        'with open("inputs.json") as f:\n'
        '    data = json.load(f)\n'
        'keyName = data["' + key_var + '"]\n'
        'cols = ' + json.dumps(col_names) + '\n'
        'numRows = len(data["' + col_names[0] + '"])\n'
        'groups = {}\n'
        'for i in range(numRows):\n'
        '    key = data[keyName][i]\n'
        '    g = groups.get(key)\n'
        '    if g is None:\n'
        '        g = {}\n'
        '        groups[key] = g\n'
        '    for col in cols:\n'
        '        val = data[col][i]\n'
        '        if isinstance(val, dict) and "path" in val:\n'
        '            val = val["path"]\n'
        '        g.setdefault(col, []).append(val)\n'
        'os.makedirs("groups", exist_ok=True)\n'
        'files = []\n'
        'for key, rows in groups.items():\n'
        '    safe = "".join(c if c.isalnum() or c == "_" else "_" for c in str(key))\n'
        '    path = f"groups/{safe}.json"\n'
        '    with open(path, "w") as f:\n'
        '        json.dump({"key": key, "rows": rows}, f)\n'
        '    files.append({"class": "File", "path": os.path.abspath(path)})\n'
        'with open("outputs.json", "w") as f:\n'
        '    json.dump(files, f)\n'
    )

    grouping_tool = {
        'id': group_tool_id,
        'class': 'CommandLineTool',
        'baseCommand': ['python3', '-c', group_script],
        'inputs': grouping_inputs,
        'outputs': [{
            'id': f'{group_tool_id}/groups',
            'type': {'type': 'array', 'items': 'File'},
            'outputBinding': {'glob': 'groups/*.json'},
        }],
        'requirements': [
            {
                'class': 'InitialWorkDirRequirement',
                'listing': listing,
            },
            {'class': 'InlineJavascriptRequirement'},
        ],
    }

    wrapper_out_names = sorted(output_schema.keys())
    wrapper_outputs = []
    for oname in wrapper_out_names:
        wrapper_outputs.append({
            'id': f'{wrapper_tool_id}/{oname}',
            'type': {'type': 'array', 'items': to_cwl_type(output_schema[oname])},
            'outputBinding': {'glob': f'{oname}s/*{oname}'},
        })

    wrapper_inputs = [{
        'id': f'{wrapper_tool_id}/group_file',
        'type': 'File',
        'inputBinding': {'position': 1},
    }]

    inner_def = step.task or {}
    dag_data = inner_def.get('dag', {})
    if dag_data and inner_def.get('class') == 'Workflow':
        original_body = _dag_to_combined_script(dag_data)
    else:
        original_body = inner_def.get('body', '')
    wrapper_script = _gen_wrapper_script(col_names, wrapper_out_names, original_body)

    wrapper_tool = {
        'id': wrapper_tool_id,
        'class': 'CommandLineTool',
        'baseCommand': ['python3', '-c', wrapper_script],
        'inputs': wrapper_inputs,
        'outputs': wrapper_outputs,
    }

    group_step_id = f'#main/group_{step_id}'
    wrap_step_id = f'#main/{step_id}'

    group_step_entry = {
        'id': group_step_id,
        'run': group_tool_id,
        'in': [{'id': f'{group_step_id}/{key_var}', 'default': group_by}],
        'out': [f'{group_step_id}/groups'],
    }
    map_source = map_info.get('source', {})
    for col in col_names:
        src = col
        if map_source.get('source') == 'table':
            cols_dict = map_source.get('columns', {})
            col_info = cols_dict.get(col, {})
            if isinstance(col_info, dict) and col_info.get('source') == 'input' and 'name' in col_info:
                src = col_info['name']
        group_step_entry['in'].append({
            'id': f'{group_step_id}/{col}',
            'source': f'#main/{src}',
        })

    wrap_step_entry = {
        'id': wrap_step_id,
        'run': wrapper_tool_id,
        'in': [{'id': f'{wrap_step_id}/group_file', 'source': f'{group_step_id}/groups'}],
        'out': [f'{wrap_step_id}/{oname}' for oname in wrapper_out_names],
        'scatter': [f'{wrap_step_id}/group_file'],
        'scatterMethod': 'dotproduct',
    }

    return [grouping_tool, wrapper_tool], [group_step_entry, wrap_step_entry]


def _gen_wrapper_script(col_names, out_names, original_body):
    import base64
    encoded_body = base64.b64encode(original_body.encode()).decode()
    lines = [
        'import base64, json, os, shutil, subprocess, sys',
        "group = json.load(open(sys.argv[1]))",
        'rows = group["rows"]',
        'cols = ' + json.dumps(col_names),
        'n = len(rows[cols[0]])',
    ]
    out_dirs = [f"'{oname}s'" for oname in out_names]
    for oname in out_names:
        lines.append(f"os.makedirs('{oname}s', exist_ok=True)")
    lines.append('for i in range(n):')
    lines.append('    row_env = os.environ.copy()')
    lines.append('    for col in cols:')
    lines.append("        row_env[col] = str(rows[col][i])")
    lines.append("    row_env['SWL_ROW_INDEX'] = str(i)")
    lines.append(f"    body = base64.b64decode('{encoded_body}').decode()")
    lines.append("    subprocess.run(body, env=row_env, shell=True, executable='/bin/bash', check=True)")
    lines.append("    out_dirs = " + json.dumps([f"{oname}s" for oname in out_names]))
    lines.append("    for f in os.listdir('.'):")
    lines.append("        if not os.path.isfile(f):")
    lines.append("            continue")
    lines.append("        if any(f.startswith(d) for d in out_dirs):")
    lines.append("            continue")
    lines.append("        for d in out_dirs:")
    lines.append("            shutil.copy(f, os.path.join(d, f))")
    lines.append("        os.remove(f)")
    return '\n'.join(lines)
