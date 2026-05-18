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
            tools.append(_tool_to_cwl(step, tool_ids[tool_id]))

    workflow = {
        'id': '#main',
        'class': 'Workflow',
        'inputs': [_workflow_input_to_cwl(workflow_id, name, spec) for name, spec in dag.inputs.items()],
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


def _step_to_cwl(workflow_id, step, tool_id):
    data = {
        'id': f'#main/{step.id}',
        'run': tool_id,
        'in': [_step_input_to_cwl(step.id, name, value) for name, value in step.bindings.items()],
        'out': [f'#main/{step.id}/{name}' for name in step.outputs],
    }
    if getattr(step, 'map', None) is not None:
        source = step.map.get('source', {})
        if source.get('source') == 'input' and 'name' in source:
            port = source['name']
            if not any(item['id'] == f'#main/{step.id}/{port}' for item in data['in']):
                data['in'].append({'id': f'#main/{step.id}/{port}', 'source': f'#main/{port}'})
            data['scatter'] = [f'#main/{step.id}/{port}']
            data['scatterMethod'] = 'dotproduct'
    return data


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
    for name, value in dag.outputs.items():
        error = _workflow_output_error(value)
        if error is not None:
            raise ValueError(f'Unsupported workflow output for CWL transpilation: {name}: {error}')

    for step in dag.steps:
        for name, value in step.bindings.items():
            error = _step_input_error(value)
            if error is not None:
                raise ValueError(f'Unsupported task input binding for CWL transpilation: {step.id}.{name}: {error}')
        for name, spec in step.task.get('outputs', {}).items():
            try:
                _interp_to_cwl_glob(spec.get('default'))
            except ValueError as exc:
                raise ValueError(f'Unsupported task output path for CWL transpilation: {step.id}.{name}: {exc}') from exc


def _step_input_error(value):
    kind = value.__class__.__name__
    if kind == 'Input':
        return None
    if kind == 'Literal':
        return None
    if kind == 'Field' and value.source.__class__.__name__ in ('TaskCall', 'StepCall', 'MappedStep'):
        return None
    if kind == 'ArrayField' and value.source.__class__.__name__ == 'MappedStep':
        return None
    if kind == 'MappedValue':
        return None
    if kind == 'Merge':
        return 'merge values are not supported'
    if kind == 'Record':
        return 'record values are not supported'
    if kind == 'Field':
        return f'field source {value.source.__class__.__name__} is not supported'
    if kind == 'ForcedFunction':
        return 'function values are not supported'
    return f'{kind} values are not supported'


def _workflow_output_error(value):
    kind = value.__class__.__name__
    if kind == 'Field' and value.source.__class__.__name__ in ('TaskCall', 'StepCall', 'MappedStep'):
        return None
    if kind == 'ArrayField' and value.source.__class__.__name__ == 'MappedStep':
        return None
    if kind == 'Literal':
        return 'literal outputs are not supported'
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
    if hasattr(value, 'name') and value.__class__.__name__ == 'Input':
        return f'#main/{value.name}'
    if value.__class__.__name__ == 'Field' and value.source.__class__.__name__ in ('TaskCall', 'StepCall', 'MappedStep'):
        return f'#main/{value.source.id}/{value.name}'
    if value.__class__.__name__ == 'ArrayField' and value.source.__class__.__name__ == 'MappedStep':
        return f'#main/{value.source.id}/{value.name}'
    if value.__class__.__name__ == 'MappedValue':
        return _mapped_value_source(value)
    if value.__class__.__name__ == 'Literal':
        return value.value
    raise ValueError(f'Unsupported binding for CWL transpilation: {value!r}')


def _infer_output_type(name, value, dag):
    if value.__class__.__name__ == 'Field' and value.source.__class__.__name__ in ('TaskCall', 'StepCall', 'MappedStep'):
        return _cwl_type(value.source.task['outputs'][name]['type'])
    if value.__class__.__name__ == 'ArrayField' and value.source.__class__.__name__ == 'MappedStep':
        return _cwl_type('[' + value.source.task['outputs'][name]['type'] + ']')
    if value.__class__.__name__ == 'MappedValue':
        return _mapped_value_type(value)
    if value.__class__.__name__ == 'Literal':
        return _cwl_type(type(value.value).__name__)
    return 'string'


def _mapped_value_source(value):
    mapped = getattr(value, 'source', None)
    element = getattr(value, 'element', None)
    if mapped is not None and getattr(mapped, '__class__', type(None)).__name__ == 'Input':
        if getattr(element, '__class__', type(None)).__name__ == 'Field':
            inner = getattr(element, 'value', None)
            if getattr(inner, '__class__', type(None)).__name__ == 'Record':
                fields = getattr(inner, 'fields', {})
                if len(fields) == 1:
                    only = next(iter(fields.values()))
                    if getattr(only, '__class__', type(None)).__name__ == 'Input':
                        return f'#main/{element.name}'
    raise ValueError(f'Unsupported mapped value for CWL transpilation: {value!r}')


def _mapped_value_type(value):
    mapped = getattr(value, 'source', None)
    if mapped is not None and getattr(mapped, '__class__', type(None)).__name__ == 'Input':
        return 'string'
    return 'string'


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
        for part in value.get('parts', []):
            if part.get('kind') == 'literal':
                parts.append(repr(part.get('text', '')))
            elif part.get('kind') == 'var':
                parts.append(f"inputs.{part.get('name')}")
            else:
                raise ValueError(f'Unsupported interpolation for CWL glob: {part!r}')
        return '$(' + ' + '.join(parts) + ')'
    raise ValueError(f'Unsupported interpolation for CWL glob: {value!r}')
