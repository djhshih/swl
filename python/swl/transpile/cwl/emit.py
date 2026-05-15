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
    for task in dag.tasks:
        tool_id = task.tool or task.id
        if tool_id not in tool_ids:
            tool_ids[tool_id] = f'#{tool_id}'
            tools.append(_tool_to_cwl(task, tool_ids[tool_id]))

    workflow = {
        'id': '#main',
        'class': 'Workflow',
        'inputs': [_workflow_input_to_cwl(workflow_id, name, spec) for name, spec in dag.inputs.items()],
        'outputs': [_workflow_output_to_cwl(workflow_id, name, value, dag) for name, value in dag.outputs.items()],
        'requirements': [],
        'steps': [_step_to_cwl(workflow_id, task, tool_ids[task.tool or task.id]) for task in dag.tasks],
    }
    return {
        'cwlVersion': 'v1.0',
        '$graph': tools + [workflow],
    }


def _tool_to_cwl(task, tool_id):
    definition = task.task or {}
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
    resource = _resource_requirement(definition.get('run', {}))
    if resource is not None:
        requirements.append(resource)
    docker = _docker_requirement(definition.get('run', {}))
    if docker is not None:
        requirements.append(docker)
    return {
        'id': tool_id,
        'class': 'CommandLineTool',
        'baseCommand': ['bash', 'script.sh'],
        'inputs': [_tool_input_to_cwl(tool_id, name, spec) for name, spec in definition.get('inputs', {}).items()],
        'outputs': [_tool_output_to_cwl(tool_id, name, spec) for name, spec in definition.get('outputs', {}).items()],
        'requirements': requirements,
    }


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


def _step_to_cwl(workflow_id, task, tool_id):
    return {
        'id': f'#main/{task.id}',
        'run': tool_id,
        'in': [_step_input_to_cwl(task.id, name, value) for name, value in task.inputs.items()],
        'out': [f'#main/{task.id}/{name}' for name in task.outputs],
    }


def _step_input_to_cwl(task_id, name, value):
    if value.__class__.__name__ == 'Literal':
        return {
            'id': f'#main/{task_id}/{name}',
            'default': value.value,
        }
    return {
        'id': f'#main/{task_id}/{name}',
        'source': _binding_source('main', value),
    }


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


def _docker_requirement(run):
    image = run.get('image', {}).get('value')
    if image is None:
        return None
    return {'class': 'DockerRequirement', 'dockerPull': image}


def _validate_supported(dag):
    for value in dag.outputs.values():
        if not _is_supported_workflow_output(value):
            raise ValueError(f'Unsupported workflow output for CWL transpilation: {value!r}')
    for task in dag.tasks:
        for value in task.inputs.values():
            if not _is_supported_step_input(value):
                raise ValueError(f'Unsupported task input binding for CWL transpilation: {value!r}')
        for name, spec in task.task.get('outputs', {}).items():
            _interp_to_cwl_glob(spec.get('default'))


def _is_supported_step_input(value):
    kind = value.__class__.__name__
    return kind in {'Input', 'Literal'} or (kind == 'Field' and value.source.__class__.__name__ == 'TaskCall')


def _is_supported_workflow_output(value):
    return value.__class__.__name__ == 'Field' and value.source.__class__.__name__ == 'TaskCall'


def _binding_source(workflow_id, value):
    if hasattr(value, 'name') and value.__class__.__name__ == 'Input':
        return f'#main/{value.name}'
    if value.__class__.__name__ == 'Field' and value.source.__class__.__name__ == 'TaskCall':
        return f'#main/{value.source.id}/{value.name}'
    if value.__class__.__name__ == 'Literal':
        return value.value
    raise ValueError(f'Unsupported binding for CWL transpilation: {value!r}')


def _infer_output_type(name, value, dag):
    if value.__class__.__name__ == 'Field' and value.source.__class__.__name__ == 'TaskCall':
        return _cwl_type(value.source.task['outputs'][name]['type'])
    if value.__class__.__name__ == 'Literal':
        return _cwl_type(type(value.value).__name__)
    return 'string'


def _cwl_type(value):
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
        for part in value.get('parts', []):
            if part.get('kind') == 'literal':
                parts.append(repr(part.get('text', '')))
            elif part.get('kind') == 'var':
                parts.append(f"inputs.{part.get('name')}")
            else:
                raise ValueError(f'Unsupported interpolation for CWL glob: {part!r}')
        return '$(' + ' + '.join(parts) + ')'
    raise ValueError(f'Unsupported interpolation for CWL glob: {value!r}')
