import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass(frozen=True)
class Input:
    name: str
    type: Optional[str] = None
    desc: Optional[str] = None


@dataclass(frozen=True)
class Literal:
    value: object


@dataclass(frozen=True)
class Record:
    fields: Dict[str, object]


@dataclass(frozen=True)
class Field:
    source: object
    name: str


@dataclass(frozen=True)
class Merge:
    left: object
    right: object


@dataclass
class TaskCall:
    id: str
    path: str
    bindings: Dict[str, object]
    outputs: List[str]
    run: Dict[str, object] = field(default_factory=dict)
    task: Optional[dict] = None
    deps: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class Output:
    name: str
    value: object


@dataclass(frozen=True)
class ForcedFunction:
    function: object
    bound: Optional[object] = None
    signature: Optional[object] = None


@dataclass
class DAG:
    inputs: Dict[str, Input]
    tasks: List[TaskCall]
    outputs: Dict[str, object]

    def to_dict(self):
        return {
            'inputs': {
                name: {
                    'type': value.type,
                    'desc': value.desc,
                }
                for name, value in self.inputs.items()
            },
            'tasks': [
                {
                    'id': task.id,
                    'path': task.path,
                    'deps': list(task.deps),
                    'inputs': dict(task.task.get('inputs', {})),
                    'bindings': {
                        name: _binding_to_binding_dict(value)
                        for name, value in task.bindings.items()
                        if not (isinstance(value, Input) and value.name == name)
                    },
                    'outputs': {
                        name: task.task['outputs'][name]
                        for name in task.outputs
                    },
                    'run': {
                        name: _run_param_to_dict(spec, task.run.get(name))
                        for name, spec in task.task.get('run', {}).items()
                    },
                    'script': task.task['body'],
                }
                for task in self.tasks
            ],
            'outputs': {name: _binding_to_dict(value) for name, value in self.outputs.items()},
        }

    @classmethod
    def from_dict(cls, data):
        inputs = {
            name: Input(name, item.get('type'), item.get('desc'))
            for name, item in data.get('inputs', {}).items()
        }
        tasks = []
        task_by_id = {}
        for item in data.get('tasks', []):
            task_outputs = item.get('outputs', {})
            task_run = item.get('run', {})
            task_inputs = item.get('inputs', {})
            task = TaskCall(
                id=item['id'],
                path=item['path'],
                bindings={},
                outputs=list(task_outputs.keys()),
                run={
                    name: _run_value_from_dict(name, spec, task_run, inputs, task_by_id)
                    for name, spec in item.get('run', {}).items()
                    if _run_value_from_dict(name, spec, task_run, inputs, task_by_id) is not None
                },
                task={
                    'doc': None,
                    'body': item.get('script', ''),
                    'inputs': task_inputs,
                    'outputs': task_outputs,
                    'run': task_run,
                },
                deps=list(item.get('deps', [])),
            )
            tasks.append(task)
            task_by_id[task.id] = task
        for task, item in zip(tasks, data.get('tasks', [])):
            task.bindings.update({
                name: inputs[name]
                for name in task.task.get('inputs', {}).keys()
                if name in inputs and name not in item.get('bindings', {})
            })
            task.bindings.update({
                name: _binding_from_binding_dict(name, value, inputs, task_by_id)
                for name, value in item.get('bindings', {}).items()
            })
        outputs = {
            name: _binding_from_dict(value, inputs, task_by_id)
            for name, value in data.get('outputs', {}).items()
        }
        return cls(inputs, tasks, outputs)

    def write(self, path):
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, sort_keys=True)

    @classmethod
    def read(cls, path):
        with open(path, 'r') as f:
            return cls.from_dict(json.load(f))


def _binding_to_dict(value):
    if isinstance(value, Input):
        return {'source': 'input', 'name': value.name}
    if isinstance(value, Literal):
        return {'source': 'literal', 'value': value.value}
    if isinstance(value, Field):
        if isinstance(value.source, TaskCall):
            return {'source': 'task', 'task': value.source.id, 'output': value.name}
        return {'source': 'field', 'field': value.name, 'value': _binding_to_dict(value.source)}
    if isinstance(value, Merge):
        return {'source': 'merge', 'left': _binding_to_dict(value.left), 'right': _binding_to_dict(value.right)}
    if isinstance(value, Record):
        return {'source': 'record', 'fields': {name: _binding_to_dict(v) for name, v in value.fields.items()}}
    if isinstance(value, TaskCall):
        return {'source': 'task_call', 'task': value.id}
    if isinstance(value, ForcedFunction):
        return {'source': 'function'}
    raise ValueError(f'Unsupported forced value for serialization: {type(value).__name__}')


def _binding_from_dict(data, inputs, tasks):
    source = data.get('source')
    if source == 'input':
        return inputs[data['name']]
    if source == 'literal':
        return Literal(data.get('value'))
    if source == 'task':
        return Field(tasks[data['task']], data['output'])
    if source == 'field':
        return Field(_binding_from_dict(data['value'], inputs, tasks), data['field'])
    if source == 'merge':
        return Merge(
            _binding_from_dict(data['left'], inputs, tasks),
            _binding_from_dict(data['right'], inputs, tasks),
        )
    if source == 'record':
        return Record({
            name: _binding_from_dict(value, inputs, tasks)
            for name, value in data.get('fields', {}).items()
        })
    if source == 'task_call':
        return tasks[data['task']]
    if source == 'function':
        return {'kind': 'function'}
    raise ValueError(f'Unsupported binding source during deserialization: {source!r}')


def _binding_to_binding_dict(value):
    if isinstance(value, Input):
        return {}
    if isinstance(value, Literal):
        return {'value': value.value}
    if isinstance(value, Field) and isinstance(value.source, TaskCall):
        return {'source': value.source.id, 'output': value.name}
    raise ValueError(f'Unsupported task binding for serialization: {type(value).__name__}')


def _binding_from_binding_dict(name, data, inputs, tasks):
    if 'value' in data:
        return Literal(data.get('value'))
    source = data.get('source')
    if source is None:
        if 'output' in data or any(key not in {'source'} for key in data.keys()):
            raise ValueError(f'Unsupported task binding during deserialization: {name}: {data!r}')
        return inputs[name]
    if 'output' not in data:
        raise ValueError(f'Unsupported task binding during deserialization: {name}: {data!r}')
    return Field(tasks[source], data['output'])


def _run_param_to_dict(spec, value):
    data = dict(spec)
    default = data.pop('default', data.pop('value', None))
    if value is None:
        data['value'] = default
        return data
    if isinstance(value, Literal):
        data['value'] = value.value
        return data
    data['value'] = _binding_to_dict(value)
    return data


def _run_value_from_dict(name, spec, defaults, inputs, tasks):
    if not isinstance(spec, dict) or 'value' not in spec:
        return None
    value = spec['value']
    default = None
    if name in defaults and isinstance(defaults[name], dict):
        default = defaults[name].get('value')
    if isinstance(value, dict) and value.get('source') is not None:
        return _binding_from_dict(value, inputs, tasks)
    if value != default:
        return Literal(value)
    return None
