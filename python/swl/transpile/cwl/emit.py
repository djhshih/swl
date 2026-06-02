import json
import os

from swl.ir.dag import DAG


def transpile_dag_file(path):
    data = json.load(open(path))
    return transpile_dag_dict(data, workflow_id='main')


def transpile_dag_dict(data, workflow_id='main'):
    dag = DAG.from_dict(data)
    _validate_supported(dag)
    tools = []
    tool_ids = {}
    for step in dag.steps:
        tool_id = step.id
        if tool_id not in tool_ids:
            tool_ids[tool_id] = f'#{tool_id}'
            tools.extend(_tool_to_cwl(step, tool_ids[tool_id]))

    workflow_inputs = dict(dag.inputs)
    for step in dag.steps:
        if getattr(step, 'map', None) is None:
            continue
        source = step.map.get('source', {})
        if source.get('source') == 'input' and 'name' in source:
            for name, typ in (getattr(step, 'input_schema', None) or {}).items():
                if name not in workflow_inputs:
                    workflow_inputs[name] = type('InputSpec', (), {'type': _as_array_type(typ), 'desc': None})()
            continue
        if source.get('source') == 'table':
            for name, typ in (getattr(step, 'input_schema', None) or {}).items():
                if name not in workflow_inputs:
                    workflow_inputs[name] = type('InputSpec', (), {'type': _as_array_type(typ), 'desc': None})()

    workflow = {
        'id': '#main',
        'class': 'Workflow',
        'inputs': [_workflow_input_to_cwl(workflow_id, name, spec) for name, spec in workflow_inputs.items()],
        'outputs': [_workflow_output_to_cwl(workflow_id, name, value, dag) for name, value in dag.outputs.items()],
        'requirements': [],
        'steps': [_step_to_cwl(workflow_id, step, tool_ids[step.id]) for step in dag.steps],
    }
    return {
        'cwlVersion': 'v1.0',
        '$graph': tools + [workflow],
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
        'type': _cwl_type(typ),
        **({'doc': desc} if desc else {}),
    }


def _workflow_output_to_cwl(workflow_id, name, value, dag):
    source = _binding_source(workflow_id, value)
    return {
        'id': f'#{workflow_id}/{name}',
        'type': _infer_output_type(name, value, dag),
        'outputSource': source,
    }


def _step_to_cwl(workflow_id, step, tool_id):
    data = {
        'id': f'#main/{step.id}',
        'run': tool_id,
        'in': [_step_input_to_cwl(step.id, name, value) for name, value in step.bindings.items()],
        'out': [f'#main/{step.id}/{name}' for name in step.outputs],
    }
    if getattr(step, 'map', None) is not None:
        ports = step.map.get('ports') or []
        if not ports:
            ports = sorted((getattr(step, 'input_schema', None) or {}).keys())
        if ports:
            source = step.map.get('source', {})
            if source.get('source') == 'input' and 'name' in source:
                for port in ports:
                    if not any(item['id'] == f'#main/{step.id}/{port}' for item in data['in']):
                        data['in'].append({'id': f'#main/{step.id}/{port}', 'source': f'#main/{port}'})
                data['scatter'] = [f'#main/{step.id}/{port}' for port in ports]
                data['scatterMethod'] = 'dotproduct'
            elif source.get('source') == 'table':
                columns = source.get('columns', {})
                for port in ports:
                    column = columns.get(port)
                    source_ref = f'#main/{port}'
                    if isinstance(column, dict) and column.get('source') == 'input' and 'name' in column:
                        source_ref = f"#main/{column['name']}"
                    if not any(item['id'] == f'#main/{step.id}/{port}' for item in data['in']):
                        data['in'].append({'id': f'#main/{step.id}/{port}', 'source': source_ref})
                data['scatter'] = [f'#main/{step.id}/{port}' for port in ports]
                data['scatterMethod'] = 'dotproduct'
    return data


def _step_input_to_cwl(task_id, name, value):
    if value.__class__.__name__ == 'Literal':
        return {
            'id': f'#main/{task_id}/{name}',
            'default': value.value,
        }
    kind, *rest = _canonical_binding(value)
    result = {
        'id': f'#main/{task_id}/{name}',
        'source': _binding_source('main', value),
    }
    if kind in ('input_field_nested', 'step_output_nested', 'tab_column_step_output_nested'):
        result['valueFrom'] = f"$({rest[-1]})"
    return result


def _tool_input_to_cwl(tool_id, name, spec):
    return {
        'id': f'{tool_id}/{name}',
        'type': _cwl_type(spec.get('type')),
        **({'doc': spec.get('desc')} if spec.get('desc') else {}),
    }


def _tool_output_to_cwl(tool_id, name, spec):
    return {
        'id': f'{tool_id}/{name}',
        'type': _cwl_type(spec.get('type')),
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


def _validate_supported(dag):
    for step in dag.steps:
        if getattr(step, 'map', None) is not None and step.map.get('group_by') is not None:
            raise ValueError(f'CWL transpilation does not yet support map_by grouping: {step.id} by {step.map.get("group_by")}')
    for name, value in dag.outputs.items():
        error = _workflow_output_error(value)
        if error is not None:
            raise ValueError(f'Unsupported workflow output for CWL transpilation: {name}: {error}')

    for step in dag.steps:
        for name, value in step.bindings.items():
            error = _step_input_error(value)
            if error is not None:
                raise ValueError(f'Unsupported step input binding for CWL transpilation: {step.id}.{name}: {error}')
        for name, spec in step.task.get('outputs', {}).items():
            try:
                _interp_to_cwl_glob(spec.get('default'))
            except ValueError as exc:
                raise ValueError(f'Unsupported step output path for CWL transpilation: {step.id}.{name}: {exc}') from exc


def _canonical_binding(value):
    kind = value.__class__.__name__
    if kind == 'Input':
        return ('input', value.name)
    if kind == 'Literal':
        return ('literal', value.value)
    if kind == 'Field' and value.source.__class__.__name__ == 'Input':
        return ('input_field', value.source.name, value.name)
    if kind == 'Field' and getattr(value.source, 'map', None) is not None:
        return ('tab_column_step_output', value.source, value.name)
    if kind == 'Field' and value.source.__class__.__name__ == 'Field':
        chain = _field_chain(value)
        root_source = chain[0].source
        root_name = chain[0].name
        field_path = '.'.join(f.name for f in chain[1:])
        if root_source.__class__.__name__ == 'Input':
            return ('input_field_nested', root_source.name, root_name, field_path)
        if getattr(root_source, 'map', None) is not None:
            return ('tab_column_step_output_nested', root_source, root_name, field_path)
        return ('step_output_nested', root_source, root_name, field_path)
    if kind == 'Field':
        return ('step_output', value.source, value.name)
    raise ValueError(f'Unsupported binding for CWL transpilation: {value!r}')


def _field_chain(value):
    chain = []
    while value.__class__.__name__ == 'Field':
        chain.append(value)
        value = value.source
    chain.reverse()
    return chain


def _step_input_error(value):
    try:
        _canonical_binding(value)
        return None
    except ValueError as exc:
        text = str(exc)
        kind = value.__class__.__name__
        if kind == 'Merge':
            return 'merge values are not supported'
        if kind == 'Record':
            return 'record values are not supported'
        if kind == 'ForcedFunction':
            return 'function values are not supported'
        if kind == 'Field':
            return f'field source {value.source.__class__.__name__} is not supported'
        return text.removeprefix('Unsupported binding for CWL transpilation: ').strip()


def _workflow_output_error(value):
    try:
        kind = _canonical_binding(value)[0]
        if kind == 'literal':
            return 'literal outputs are not supported'
        return None
    except ValueError:
        kind = value.__class__.__name__
        if kind == 'Merge':
            return 'merge outputs are not supported'
        if kind == 'Record':
            return 'record outputs are not supported'
        if kind == 'ForcedFunction':
            return 'function outputs are not supported'
        if kind == 'Field':
            return f'field source {value.source.__class__.__name__} is not supported'
        return f'{kind} outputs are not supported'


def _binding_source(workflow_id, value):
    kind, *rest = _canonical_binding(value)
    if kind == 'input':
        return f'#main/{rest[0]}'
    if kind == 'step_output':
        step, output = rest
        return f'#main/{step.id}/{output}'
    if kind == 'tab_column_step_output':
        step, output = rest
        return f'#main/{step.id}/{output}'
    if kind == 'input_field':
        input_name, field_name = rest
        return f'#main/{input_name}/{field_name}'
    if kind == 'input_field_nested':
        input_name, root_name, field_path = rest
        return f'#main/{input_name}'
    if kind == 'step_output_nested':
        step, root_name, field_path = rest
        return f'#main/{step.id}/{root_name}'
    if kind == 'tab_column_step_output_nested':
        step, root_name, field_path = rest
        return f'#main/{step.id}/{root_name}'
    if kind == 'literal':
        return rest[0]
    raise ValueError(f'Unsupported binding for CWL transpilation: {value!r}')


def _infer_output_type(name, value, dag):
    kind, *rest = _canonical_binding(value)
    if kind == 'step_output':
        step, output = rest
        return _cwl_type(step.task['outputs'][output]['type'])
    if kind == 'tab_column_step_output':
        step, output = rest
        return _cwl_type('[' + step.task['outputs'][output]['type'] + ']')
    if kind in ('step_output_nested', 'tab_column_step_output_nested'):
        step, root_name, field_path = rest
        return _cwl_type(step.task['outputs'][root_name]['type'])
    if kind in ('input_field', 'input_field_nested'):
        return 'string'
    if kind == 'literal':
        return _cwl_type(type(rest[0]).__name__)
    return 'string'



def _as_array_type(value):
    if value is None or (isinstance(value, str) and value.startswith('[') and value.endswith(']')):
        return value
    return f'[{value}]'


def _cwl_type(value):
    if value == '[file]':
        return {'type': 'array', 'items': 'File'}
    if value == '[str]':
        return {'type': 'array', 'items': 'string'}
    if value == '[int]':
        return {'type': 'array', 'items': 'int'}
    if value == '[float]':
        return {'type': 'array', 'items': 'float'}
    return {
        'file': 'File',
        'str': 'string',
        'string': 'string',
        'int': 'int',
        'float': 'float',
        'bool': 'boolean',
        'boolean': 'boolean',
    }.get(value, 'string')


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
