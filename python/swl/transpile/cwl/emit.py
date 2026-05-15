import json
import os

from swl.ir.dag import DAG


def transpile_dag_file(path):
    data = json.load(open(path))
    workflow_id = os.path.splitext(os.path.basename(path))[0]
    return transpile_dag_dict(data, workflow_id=workflow_id)


def transpile_dag_dict(data, workflow_id='main'):
    dag = DAG.from_dict(data)
    tools = []
    tool_ids = {}
    for task in dag.tasks:
        tool_id = task.name
        if tool_id not in tool_ids:
            tool_ids[tool_id] = f'#{tool_id}'
            tools.append(_tool_to_cwl(task, tool_ids[tool_id]))

    workflow = {
        'id': f'#{workflow_id}',
        'class': 'Workflow',
        'inputs': [_workflow_input_to_cwl(workflow_id, name, spec) for name, spec in dag.inputs.items()],
        'outputs': [_workflow_output_to_cwl(workflow_id, name, value, dag) for name, value in dag.outputs.items()],
        'requirements': [],
        'steps': [_step_to_cwl(workflow_id, task, tool_ids[task.name]) for task in dag.tasks],
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
        'id': f'#{workflow_id}/{task.id}',
        'run': tool_id,
        'in': [
            {
                'id': f'#{workflow_id}/{task.id}/{name}',
                'source': _binding_source(workflow_id, value),
            }
            for name, value in task.inputs.items()
        ],
        'out': [f'#{workflow_id}/{task.id}/{name}' for name in task.outputs],
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


def _binding_source(workflow_id, value):
    if hasattr(value, 'name') and value.__class__.__name__ == 'Input':
        return f'#{workflow_id}/{value.name}'
    if value.__class__.__name__ == 'Field' and value.source.__class__.__name__ == 'TaskCall':
        return f'#{workflow_id}/{value.source.id}/{value.name}'
    raise ValueError(f'Unsupported binding for CWL transpilation: {value!r}')


def _infer_output_type(name, value, dag):
    if value.__class__.__name__ == 'Field' and value.source.__class__.__name__ == 'TaskCall':
        return _cwl_type(value.source.task['outputs'][name]['type'])
    return 'string'


def _cwl_type(value):
    return {
        'file': 'File',
        'str': 'string',
        'int': 'int',
        'float': 'float',
        'bool': 'boolean',
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
