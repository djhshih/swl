import json
import os

from swl.dag.node import DAG, Field, Input, Literal, OutputSpec, Record, StepCall
from swl.types import to_array_type, to_cwl_type


def transpile_dag_file(path):
    data = json.load(open(path))
    return transpile_dag_dict(data, workflow_id='main')


def transpile_dag_dict(data, workflow_id='main'):
    dag = DAG.from_dict(data)
    dag.validate()
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
                tool, step_entry, source = _emit_record_tool(step.id, name, value, dag)
                record_tools.append(tool)
                record_map[(step.id, name)] = (step_entry, source)

    for name, output in list(dag.outputs.items()):
        value = output.value if isinstance(output, OutputSpec) else output
        if isinstance(value, Record):
            tool, step_entry, source = _emit_record_tool('outputs', name, value, dag)
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
        for name, typ in (getattr(step, 'input_schema', None) or {}).items():
            if name not in workflow_inputs:
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
            outputs.append(_workflow_output_to_cwl(workflow_id, name, output))

    workflow = {
        'id': '#main',
        'class': 'Workflow',
        'inputs': [_workflow_input_to_cwl(workflow_id, name, spec) for name, spec in workflow_inputs.items()],
        'outputs': outputs,
        'requirements': [{'class': 'ScatterFeatureRequirement'}, {'class': 'MultipleInputFeatureRequirement'}],
        'steps': wf_steps,
    }
    return {
        'cwlVersion': 'v1.0',
        '$graph': all_tools + [workflow],
    }


def _tool_to_cwl(step, tool_id):
    definition = step.task or {}
    if definition.get('class') == 'Workflow':
        packed = transpile_dag_dict(definition['dag'], workflow_id=tool_id[1:])
        graph = []
        for item in packed.get('$graph', []):
            copied = dict(item)
            if copied.get('id') == '#main':
                copied['id'] = tool_id
            elif copied.get('id', '').startswith('#main/'):
                copied['id'] = copied['id'].replace('#main/', f'{tool_id}/', 1)
            graph.append(copied)
        return graph
    run = definition.get('run', {})
    requirements = [
        {
            'class': 'InitialWorkDirRequirement',
            'listing': [
                {
                    'entryname': 'script.sh',
                    'entry': definition.get('body', ''),
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
    if any(_has_expr_interpolation(spec) for spec in definition.get('outputs', {}).values()):
        requirements.append({'class': 'InlineJavascriptRequirement'})
    tool = {
        'id': tool_id,
        'class': 'CommandLineTool',
        'baseCommand': ['bash', 'script.sh'],
        'inputs': [_tool_input_to_cwl(tool_id, name, spec) for name, spec in definition.get('inputs', {}).items()],
        'outputs': [_tool_output_to_cwl(tool_id, name, spec) for name, spec in definition.get('outputs', {}).items()],
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


def _workflow_output_to_cwl(workflow_id, name, output):
    source = _binding_source(output.value)
    return {
        'id': f'#{workflow_id}/{name}',
        'type': to_cwl_type(output.type),
        'outputSource': source,
    }


def _step_to_cwl(workflow_id, step, tool_id, record_map=None):
    if record_map is None:
        record_map = {}
    inputs = []
    for name, value in step.bindings.items():
        if (step.id, name) in record_map:
            _, source = record_map[(step.id, name)]
            inputs.append({
                'id': f'#main/{step.id}/{name}',
                'source': source,
            })
        else:
            inputs.append(_step_input_to_cwl(step.id, name, value))
    data = {
        'id': f'#main/{step.id}',
        'run': tool_id,
        'in': inputs,
        'out': [f'#main/{step.id}/{name}' for name in step.outputs],
    }
    if getattr(step, 'map', None) is not None:
        scatter_ports = step.map.get('scatter') or sorted((getattr(step, 'input_schema', None) or {}).keys())
        if scatter_ports:
            source = step.map.get('source', {})
            if source.get('source') == 'input' and 'name' in source:
                for port in scatter_ports:
                    if not any(item['id'] == f'#main/{step.id}/{port}' for item in data['in']):
                        data['in'].append({'id': f'#main/{step.id}/{port}', 'source': f'#main/{port}'})
            elif source.get('source') == 'table':
                columns = source.get('columns', {})
                for port in scatter_ports:
                    column = columns.get(port)
                    source_ref = f'#main/{port}'
                    if isinstance(column, dict) and column.get('source') == 'input' and 'name' in column:
                        source_ref = f"#main/{column['name']}"
                    if not any(item['id'] == f'#main/{step.id}/{port}' for item in data['in']):
                        data['in'].append({'id': f'#main/{step.id}/{port}', 'source': source_ref})
            data['scatter'] = [f'#main/{step.id}/{port}' for port in scatter_ports]
            data['scatterMethod'] = 'dotproduct'
    return data


def _step_input_to_cwl(task_id, name, value):
    if isinstance(value, Literal):
        return {
            'id': f'#main/{task_id}/{name}',
            'default': value.value,
        }
    result = {
        'id': f'#main/{task_id}/{name}',
        'source': _binding_source(value),
    }
    field_path = _field_path_after_first(value)
    if field_path is not None:
        result['valueFrom'] = f"$({field_path})"
    return result


def _tool_input_to_cwl(tool_id, name, spec):
    return {
        'id': f'{tool_id}/{name}',
        'type': to_cwl_type(spec.get('type')),
        **({'doc': spec.get('desc')} if spec.get('desc') else {}),
    }


def _tool_output_to_cwl(tool_id, name, spec):
    return {
        'id': f'{tool_id}/{name}',
        'type': to_cwl_type(spec.get('type')),
        'outputBinding': {
            'glob': _interp_to_cwl_glob(spec.get('default')),
        },
        **({'doc': spec.get('desc')} if spec.get('desc') else {}),
    }


def _resource_requirement(run):
    req = {'class': 'ResourceRequirement'}
    if 'cpu' in run:
        req['coresMin'] = run['cpu'].get('value')
    if 'memory' in run:
        req['ramMin'] = run['memory'].get('value')
    if len(req) == 1:
        return None
    return req


def _hints_from_run(run):
    hints = []
    if 'time' in run:
        hints.append({'class': 'TimeLimit', 'timeLimit': run['time'].get('value')})
    return hints


def _docker_requirement(run):
    image = run.get('image', {}).get('value')
    if image is None:
        return None
    return {'class': 'DockerRequirement', 'dockerPull': image}


def _field_path_after_first(value):
    if not isinstance(value, Field):
        return None
    chain = []
    v = value
    while isinstance(v, Field):
        chain.append(v)
        v = v.source
    if len(chain) <= 1:
        return None
    chain.reverse()
    return '.'.join(f.name for f in chain[1:])


def _binding_source(value):
    if isinstance(value, Input):
        return f'#main/{value.name}'
    if isinstance(value, Literal):
        return value.value
    if isinstance(value, StepCall):
        return f'#main/{value.id}'
    if isinstance(value, Field):
        src = value.source
        if isinstance(src, Input):
            return f'#main/{src.name}/{value.name}'
        if isinstance(src, StepCall):
            return f'#main/{src.id}/{value.name}'
        if isinstance(src, Field):
            root, root_name, _ = _field_chain_parts(value)
            if isinstance(root, Input):
                return f'#main/{root.name}'
            if isinstance(root, StepCall):
                return f'#main/{root.id}/{root_name}'
    raise ValueError(f'Unsupported binding for CWL transpilation: {value!r}')


def _field_chain_parts(value):
    chain = []
    v = value
    while isinstance(v, Field):
        chain.append(v)
        v = v.source
    if not chain:
        return value, None, None
    chain.reverse()
    return chain[0].source, chain[0].name, '.'.join(f.name for f in chain[1:])


def _validate_supported(dag):
    for step in dag.steps:
        if getattr(step, 'map', None) is not None and step.map.get('group_by') is not None:
            _validate_map_by_preconditions(step)
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        if isinstance(value, Literal):
            raise ValueError(f'Unsupported workflow output for CWL transpilation: {name}: literal outputs are not supported')
    for step in dag.steps:
        for name, spec in step.task.get('outputs', {}).items():
            try:
                _interp_to_cwl_glob(spec.get('default'))
            except ValueError as exc:
                raise ValueError(f'Unsupported step output path for CWL transpilation: {step.id}.{name}: {exc}') from exc


def _interp_to_cwl_glob(value):
    if value is None:
        return '*'
    if value.get('kind') == 'word':
        parts = []
        has_expr = False
        for part in value.get('parts', []):
            if part.get('kind') == 'literal':
                parts.append(repr(part.get('text', '')))
            elif part.get('kind') == 'var':
                parts.append(f"inputs.{part.get('name')}")
            elif part.get('kind') == 'expr':
                parts.append(f"({part.get('text')})")
                has_expr = True
            else:
                raise ValueError(f'Unsupported interpolation for CWL glob: {part!r}')
        expr = '$(' + ' + '.join(parts) + ')'
        return expr
    raise ValueError(f'Unsupported interpolation for CWL glob: {value!r}')


def _has_expr_interpolation(spec):
    default = spec.get('default') if isinstance(spec, dict) else getattr(spec, 'default', None)
    if default is None:
        return False
    if isinstance(default, dict) and default.get('kind') == 'word':
        return any(part.get('kind') == 'expr' for part in default.get('parts', []))
    return False


# record bindings ---------------------------------------------------------

def _infer_record_field_type(field_value, dag):
    if isinstance(field_value, Input):
        spec = dag.inputs.get(field_value.name)
        typ = spec.type if hasattr(spec, 'type') else (spec or {}).get('type') if isinstance(spec, dict) else None
        return to_cwl_type(typ)
    if isinstance(field_value, Literal):
        return to_cwl_type(type(field_value.value).__name__)
    if isinstance(field_value, StepCall):
        return 'string'
    if isinstance(field_value, Field):
        root, root_name, _ = _field_chain_parts(field_value)
        if isinstance(root, StepCall) and root.task:
            out_spec = root.task.get('outputs', {}).get(root_name, {})
            if out_spec:
                return to_cwl_type(out_spec.get('type') or 'str')
        return 'string'
    return 'string'


def _emit_record_tool(step_id, binding_name, record, dag):
    tool_id = f'#rec_{step_id}_{binding_name}'
    step_entry_id = f'#main/rec_{step_id}_{binding_name}'
    field_names = sorted(record.fields.keys())
    tool_inputs = []
    step_inputs = []
    for fname in field_names:
        fvalue = record.fields[fname]
        typ = _infer_record_field_type(fvalue, dag)
        tool_inputs.append({'id': f'{tool_id}/{fname}', 'type': typ})
        if isinstance(fvalue, Literal):
            step_inputs.append({'id': f'{step_entry_id}/{fname}', 'default': fvalue.value})
        elif isinstance(fvalue, (Input, Field, StepCall)):
            step_inputs.append({'id': f'{step_entry_id}/{fname}', 'source': _binding_source(fvalue)})
        else:
            step_inputs.append({'id': f'{step_entry_id}/{fname}', 'default': None})
    js_lines = ['var fs = require("fs");']
    for fname in field_names:
        js_lines.append(f'record.{fname} = inputs.{fname};')
    js_lines.append('fs.writeFileSync("record.json", JSON.stringify(record));')
    js_lines.append('return {"class": "File", "path": "record.json"};')
    js_expression = 'var record = {};\n' + '\n'.join(js_lines)
    tool = {
        'id': tool_id,
        'class': 'ExpressionTool',
        'requirements': [{'class': 'InlineJavascriptRequirement'}],
        'inputs': tool_inputs,
        'outputs': [{
            'id': f'{tool_id}/record_out',
            'type': 'File',
            'outputBinding': {'glob': 'record.json'},
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


def _emit_map_by_graph(step, dag):
    map_info = step.map or {}
    group_by = map_info.get('group_by')
    schema = step.input_schema or {}
    output_schema = step.output_schema or {}
    step_id = step.id
    group_tool_id = f'#group_{step_id}'
    wrapper_tool_id = f'#wrap_{step_id}'

    col_names = sorted(schema.keys())
    js_inputs = {}
    for col in col_names:
        js_inputs[col] = f"inputs.{col}[i]"
    key_var = 'group_by_key'

    group_js = (
        'var fs = require("fs");\n'
        'var keyName = inputs.' + key_var + ';\n'
        'var cols = ' + json.dumps(col_names) + ';\n'
        'var numRows = inputs.' + col_names[0] + '.length;\n'
        'var groups = {};\n'
        'for (var i = 0; i < numRows; i++) {\n'
        '    var key = inputs[keyName][i];\n'
        '    var g = groups[key];\n'
        '    if (!g) { g = {}; groups[key] = g; }\n'
        '    for (var ci = 0; ci < cols.length; ci++) {\n'
        '        var col = cols[ci];\n'
        '        if (!g[col]) g[col] = [];\n'
        '        g[col].push(inputs[col][i]);\n'
        '    }\n'
        '}\n'
        'try { fs.mkdirSync("groups"); } catch(e) {}\n'
        'var files = [];\n'
        'for (var key in groups) {\n'
        '    var safe = key.replace(/[^a-zA-Z0-9_]/g, "_");\n'
        '    var path = "groups/" + safe + ".json";\n'
        '    fs.writeFileSync(path, JSON.stringify({key: key, rows: groups[key]}));\n'
        '    files.push({"class": "File", "path": path});\n'
        '}\n'
        'return files;\n'
    )

    grouping_tool = {
        'id': group_tool_id,
        'class': 'ExpressionTool',
        'requirements': [{'class': 'InlineJavascriptRequirement'}],
        'inputs': [{  # group_by_key input
            'id': f'{group_tool_id}/{key_var}',
            'type': 'string',
        }],
        'outputs': [{
            'id': f'{group_tool_id}/groups',
            'type': {'type': 'array', 'items': 'File'},
            'outputBinding': {'glob': 'groups/*.json'},
        }],
        'expression': group_js,
    }
    for col in col_names:
        grouping_tool['inputs'].append({
            'id': f'{group_tool_id}/{col}',
            'type': {'type': 'array', 'items': to_cwl_type(schema[col])},
        })

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
    }]

    inner_def = step.task or {}
    original_body = inner_def.get('body', '')
    row_loop = _gen_wrapper_script(col_names, wrapper_out_names, original_body)

    wrapper_tool = {
        'id': wrapper_tool_id,
        'class': 'CommandLineTool',
        'baseCommand': ['bash', 'wrapper.sh'],
        'inputs': wrapper_inputs,
        'outputs': wrapper_outputs,
        'requirements': [
            {
                'class': 'InitialWorkDirRequirement',
                'listing': [{'entryname': 'wrapper.sh', 'entry': row_loop}],
            },
            {'class': 'InlineJavascriptRequirement'},
        ],
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
    lines = ['#!/bin/bash', 'set -e', 'group_file="$1"', '']
    lines.append('python3 -c "')
    lines.append('import json, os, subprocess, sys')
    lines.append("group = json.load(open('${group_file}'))")
    lines.append('rows = group[\"rows\"]')
    lines.append('cols = ' + json.dumps(col_names))
    lines.append('n = len(rows[cols[0]])')
    for oname in out_names:
        lines.append(f"os.makedirs('{oname}s', exist_ok=True)")
    lines.append('for i in range(n):')
    lines.append('    row_env = os.environ.copy()')
    lines.append('    for col in cols:')
    lines.append('        row_env[col] = str(rows[col][i])')
    lines.append("    row_env['SWL_ROW_INDEX'] = str(i)")
    lines.append("    subprocess.run(['bash', '-c', '''")
    for line in original_body.split('\n'):
        lines.append('    ' + line)
    lines.append("    '''], env=row_env, shell=True, check=True)")
    for oname in out_names:
        lines.append(f"    for f in os.listdir('.'):")
        lines.append(f"        if os.path.isfile(f) and not f.startswith('{oname}s'):")
        lines.append(f"            os.rename(f, '{oname}s/' + f)")
    lines.append('"')
    return '\n'.join(lines)
