import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

"""
Interpolation JSON schema for task parameter defaults.

A default value is serialized as one of these forms:

Word (interpolation with parts):
    {"kind": "word", "parts": [
        {"kind": "literal", "text": "path/to/"},
        {"kind": "var", "name": "sample_id"},
        {"kind": "literal", "text": ".bam"}
    ]}

Literal (plain text):
    {"kind": "literal", "text": "42"}

Var (variable reference — resolved from executor bindings):
    {"kind": "var", "name": "outbase"}

Expr (arbitrary expression — executor may shell-evaluate or reject):
    {"kind": "expr", "text": "sample_id + \"_out\""}

Resolution rules:
- Literal: used as-is.
- Var: resolved from the executor's runtime bindings dict by name.
- Expr: the executor MAY evaluate the text as a shell expression, or reject
  it if evaluation is not supported. The compiler does not validate Expr
  contents beyond parsing.
- Word: concatenate the resolved parts in order.
"""


@dataclass(frozen=True)
class Input:
    name: str
    type: Optional[str] = None
    desc: Optional[str] = None
    optional: bool = False


@dataclass(frozen=True)
class Literal:
    value: object


@dataclass(frozen=True)
class Record:
    fields: Dict[str, object]


@dataclass(frozen=True)
class TableSource:
    name: str
    columns: Dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class Field:
    source: object
    name: str


@dataclass(frozen=True)
class Merge:
    left: object
    right: object


@dataclass
class StepCall:
    id: str
    path: str
    bindings: Dict[str, object]
    outputs: List[str]
    source: object = None
    run: Dict[str, object] = field(default_factory=dict)
    task: Optional[dict] = None
    deps: List[str] = field(default_factory=list)
    type: str = 'task'
    map: Optional[dict] = None
    input_schema: Optional[dict] = None
    output_schema: Optional[dict] = None


@dataclass(frozen=True)
class Output:
    name: str
    value: object


@dataclass(frozen=True)
class ForcedFunction:
    function: object
    bound: Optional[object] = None
    signature: Optional[object] = None
    satisfied: Set[str] = field(default_factory=set)


@dataclass(init=False)
class DAG:
    inputs: Dict[str, Input]
    steps: List[StepCall]
    outputs: Dict[str, object]

    def __init__(self, inputs, steps=None, outputs=None):
        self.inputs = inputs
        self.steps = steps if steps is not None else []
        self.outputs = outputs if outputs is not None else {}

    def to_dict(self):
        return {
            'inputs': {
                name: {
                    'type': value.type,
                    'desc': value.desc,
                    'optional': value.optional,
                }
                for name, value in self.inputs.items()
            },
            'steps': [
                {
                    'id': step.id,
                    'type': step.type,
                    'path': step.path,
                    **({'map': step.map} if step.map is not None else {}),
                    **({'input_schema': step.input_schema} if step.input_schema is not None else {}),
                    **({'output_schema': step.output_schema} if step.output_schema is not None else {}),
                    'deps': list(step.deps),
                    'inputs': dict(step.task.get('inputs', {})),
                    'bindings': {
                        name: _binding_to_binding_dict(value)
                        for name, value in step.bindings.items()
                        if not (isinstance(value, Input) and value.name == name)
                    },
                    'outputs': {
                        name: step.task['outputs'][name]
                        for name in step.outputs
                    },
                    'run': {
                        name: _run_param_to_dict(spec, step.run.get(name))
                        for name, spec in step.task.get('run', {}).items()
                    },
                    'script': step.task['body'],
                    **({'definition': step.task} if step.task.get('class') == 'Workflow' else {}),
                }
                for step in self.steps
            ],
            'outputs': {name: _binding_to_dict(value) for name, value in self.outputs.items()},
        }

    def validate(self):
        step_ids = {step.id for step in self.steps}
        for step in self.steps:
            for dep in step.deps:
                if dep not in step_ids:
                    raise ValueError(f'Step {step.id} depends on unknown step: {dep}')
        visited = set()
        in_stack = set()
        def visit(step_id):
            if step_id in in_stack:
                raise ValueError(f'Circular dependency detected in DAG involving step: {step_id}')
            if step_id in visited:
                return
            visited.add(step_id)
            in_stack.add(step_id)
            step = next(s for s in self.steps if s.id == step_id)
            for dep in step.deps:
                visit(dep)
            in_stack.remove(step_id)
        for step in self.steps:
            if step.id not in visited:
                visit(step.id)

    @classmethod
    def from_dict(cls, data):
        inputs = {
            name: Input(name, item.get('type'), item.get('desc'), item.get('optional', False))
            for name, item in data.get('inputs', {}).items()
        }
        steps_data = data.get('steps', [])
        steps = []
        step_by_id = {}
        for item in steps_data:
            step_outputs = item.get('outputs', {})
            step_run = item.get('run', {})
            step_inputs = item.get('inputs', {})
            step_cls = StepCall
            step = step_cls(
                id=item['id'],
                path=item['path'],
                bindings={},
                outputs=list(step_outputs.keys()),
                run={
                    name: value
                    for name, spec in item.get('run', {}).items()
                    if (value := _run_value_from_dict(name, spec, step_run, inputs, step_by_id)) is not None
                },
                task=item.get('definition', {
                    'doc': None,
                    'body': item.get('script', ''),
                    'inputs': step_inputs,
                    'outputs': step_outputs,
                    'run': step_run,
                }),
                deps=list(item.get('deps', [])),
                type=item.get('type', 'task'),
                **({'source': _binding_from_dict(item['map']['source'], inputs, step_by_id) if item.get('map', {}).get('source') is not None or item.get('map', {}).get('fields') is not None else None, 'map': item.get('map'), 'input_schema': item.get('input_schema'), 'output_schema': item.get('output_schema')} if item.get('map') is not None else {}),
            )
            steps.append(step)
            step_by_id[step.id] = step
        for step, item in zip(steps, steps_data):
            step.bindings.update({
                name: inputs[name]
                for name in step.task.get('inputs', {}).keys()
                if name in inputs and name not in item.get('bindings', {})
            })
            step.bindings.update({
                name: _binding_from_binding_dict(name, value, inputs, step_by_id)
                for name, value in item.get('bindings', {}).items()
            })
        outputs = {
            name: _binding_from_dict(value, inputs, step_by_id)
            for name, value in data.get('outputs', {}).items()
        }
        return cls(inputs, steps, outputs)

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
        if isinstance(value.source, StepCall):
            return {'step': value.source.id, 'output': value.name}
        return {'source': 'field', 'field': value.name, 'value': _binding_to_dict(value.source)}
    if isinstance(value, Merge):
        return {'source': 'merge', 'left': _binding_to_dict(value.left), 'right': _binding_to_dict(value.right)}
    if isinstance(value, Record):
        return {'source': 'record', 'fields': {name: _binding_to_dict(v) for name, v in value.fields.items()}}
    if isinstance(value, TableSource):
        return {'source': 'table', 'name': value.name, 'columns': {name: _binding_to_dict(v) for name, v in value.columns.items()}}
    if isinstance(value, StepCall):
        return {'source': 'step_call', 'step': value.id}
    if isinstance(value, ForcedFunction):
        return {'source': 'function'}
    raise ValueError(f'Unsupported forced value for serialization: {type(value).__name__}')


def _binding_from_dict(data, inputs, steps):
    if 'step' in data and 'output' in data:
        step = steps[data['step']]
        return Field(step, data['output'])
    source = data.get('source')
    if source == 'input':
        name = data['name']
        return inputs.get(name, Input(name, optional=False))
    if source == 'value':
        return _binding_from_dict(data['value'], inputs, steps)
    if source == 'literal':
        return Literal(data.get('value'))
    if source == 'field':
        return Field(_binding_from_dict(data['value'], inputs, steps), data['field'])
    if source == 'merge':
        return Merge(
            _binding_from_dict(data['left'], inputs, steps),
            _binding_from_dict(data['right'], inputs, steps),
        )
    if source == 'record':
        return Record({
            name: _binding_from_dict(value, inputs, steps)
            for name, value in data.get('fields', {}).items()
        })
    if source == 'table':
        return TableSource(
            data['name'],
            {
                name: _binding_from_dict(value, inputs, steps)
                for name, value in data.get('columns', {}).items()
            },
        )
    if source == 'step':
        step = steps[data['step']]
        return Field(step, data['output'])
    if source == 'step_call':
        return steps[data['step']]
    if source == 'function':
        return {'kind': 'function'}
    raise ValueError(f'Unsupported binding source during deserialization: {source!r}')


def _binding_to_binding_dict(value):
    if isinstance(value, Input):
        return {}
    if isinstance(value, Literal):
        return {'value': value.value}
    if isinstance(value, Record):
        return {'source': 'record', 'fields': {name: _binding_to_dict(v) for name, v in value.fields.items()}}
    if isinstance(value, Field) and isinstance(value.source, StepCall):
        return {'source': value.source.id, 'output': value.name}
    if isinstance(value, Field) and isinstance(value.source, Input):
        return {'kind': 'field', 'source': {'source': 'input', 'name': value.source.name}, 'field': value.name}
    if isinstance(value, Field) and isinstance(value.source, Field):
        return _nested_field_to_dict(value)
    raise ValueError(f'Unsupported step binding for serialization: {type(value).__name__}')


def _nested_field_to_dict(field):
    chain = []
    while isinstance(field, Field) and isinstance(field.source, Field):
        chain.append(field)
        field = field.source
    inner = _binding_to_binding_dict(field)
    if isinstance(inner, dict) and 'source' in inner and 'output' in inner:
        inner = {'step': inner['source'], 'output': inner['output']}
    result = inner
    for f in reversed(chain):
        result = {'kind': 'field', 'source': result, 'field': f.name}
    return result


def _binding_from_binding_dict(name, data, inputs, steps):
    if 'value' in data:
        return Literal(data.get('value'))
    if data.get('kind') == 'field':
        return Field(_binding_from_dict(data['source'], inputs, steps), data['field'])
    if data.get('source') == 'record':
        return Record({
            fname: _binding_from_dict(fdata, inputs, steps)
            for fname, fdata in data.get('fields', {}).items()
        })
    source = data.get('source')
    if source is None:
        if 'output' in data or any(key not in {'source'} for key in data.keys()):
            raise ValueError(f'Unsupported step binding during deserialization: {name}: {data!r}')
        return inputs[name]
    if 'output' not in data:
        raise ValueError(f'Unsupported step binding during deserialization: {name}: {data!r}')
    step = steps[source]
    return Field(step, data['output'])


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


def _run_value_from_dict(name, spec, defaults, inputs, steps):
    if not isinstance(spec, dict) or 'value' not in spec:
        return None
    value = spec['value']
    default = None
    if name in defaults and isinstance(defaults[name], dict):
        default = defaults[name].get('value')
    if isinstance(value, dict) and value.get('source') is not None:
        return _binding_from_dict(value, inputs, steps)
    if value != default:
        return Literal(value)
    return None


def _resolve_interp_part(part, bindings):
    kind = part.get('kind')
    if kind == 'literal':
        return part['text']
    if kind == 'var':
        name = part['name']
        if name in bindings:
            return str(bindings[name])
        return '${' + name + '}'
    if kind == 'expr':
        return '${' + part['text'] + '}'
    return ''


def resolve_default(default, bindings=None):
    if bindings is None:
        bindings = {}
    if not isinstance(default, dict):
        return default
    kind = default.get('kind')
    if kind == 'literal':
        return default['text']
    if kind == 'var':
        name = default['name']
        if name in bindings:
            return str(bindings[name])
        return '${' + name + '}'
    if kind == 'expr':
        return '${' + default['text'] + '}'
    if kind == 'word':
        parts = default.get('parts', [])
        return ''.join(_resolve_interp_part(part, bindings) for part in parts)
    return None
