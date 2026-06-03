import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

from swl.dag.binding import binding_from_dict, binding_to_dict

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

Var (variable reference -- resolved from executor bindings):
    {"kind": "var", "name": "outbase"}

Expr (arbitrary expression -- executor may shell-evaluate or reject):
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
class OutputSpec:
    type: Optional[str]
    value: object
    desc: Optional[str] = None
    optional: bool = False


@dataclass(frozen=True)
class ForcedFunction:
    function: object
    bound: Optional[object] = None
    signature: Optional[object] = None
    satisfied: Set[str] = field(default_factory=set)


@dataclass
class DAG:
    inputs: Dict[str, Input]
    steps: List[StepCall] = field(default_factory=list)
    outputs: Dict[str, OutputSpec] = field(default_factory=dict)

    def to_dict(self):
        return {
            'inputs': {
                name: {
                    'type': value.type,
                    'desc': value.desc,
                    **({'optional': True} if value.optional else {}),
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
                        name: binding_to_dict(value)
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
            'outputs': {
                name: (
                    {
                        'type': value.type,
                        'desc': value.desc,
                        **({'optional': True} if value.optional else {}),
                        'value': binding_to_dict(value.value),
                    }
                    if isinstance(value, OutputSpec) or (hasattr(value, 'value') and hasattr(value, 'type'))
                    else binding_to_dict(value)
                )
                for name, value in self.outputs.items()
            },
        }

    def validate(self):
        step_by_id = {step.id: step for step in self.steps}
        self._validate_step_deps(step_by_id)
        self._validate_acyclic(step_by_id)
        for step in self.steps:
            if step.map is not None:
                self._validate_mapped_ports(step)
        self._check_no_disallowed_values()
        self._check_output_types()

    def _validate_step_deps(self, step_by_id):
        for step in self.steps:
            for dep in step.deps:
                if dep not in step_by_id:
                    raise ValueError(f'Step {step.id} depends on unknown step: {dep}')

    def _validate_acyclic(self, step_by_id):
        visited = set()
        in_stack = set()

        def visit(step_id):
            if step_id in in_stack:
                raise ValueError(f'Circular dependency detected in DAG involving step: {step_id}')
            if step_id in visited:
                return
            visited.add(step_id)
            in_stack.add(step_id)
            for dep in step_by_id[step_id].deps:
                visit(dep)
            in_stack.remove(step_id)

        for step_id in step_by_id:
            visit(step_id)

    @staticmethod
    def _walk_values(value):
        if isinstance(value, Merge):
            yield from DAG._walk_values(value.left)
            yield from DAG._walk_values(value.right)
        elif isinstance(value, Record):
            for fv in value.fields.values():
                yield from DAG._walk_values(fv)
        elif isinstance(value, Field):
            yield from DAG._walk_values(value.source)
        elif isinstance(value, ForcedFunction):
            if value.bound is not None:
                yield from DAG._walk_values(value.bound)

    @staticmethod
    def _check_value(value, label):
        if isinstance(value, Merge):
            raise ValueError(f'{label} contains a Merge value')
        if isinstance(value, ForcedFunction):
            raise ValueError(f'{label} contains a ForcedFunction')
        for v in DAG._walk_values(value):
            if isinstance(v, Merge):
                raise ValueError(f'{label} contains a Merge value nested inside {type(v).__name__}')
            if isinstance(v, ForcedFunction):
                raise ValueError(f'{label} contains a ForcedFunction nested inside {type(v).__name__}')

    def _check_no_disallowed_values(self):
        for step in self.steps:
            for name, value in step.bindings.items():
                DAG._check_value(value, f'Step {step.id} binding {name!r}')
        for name, output in self.outputs.items():
            value = output.value if isinstance(output, OutputSpec) else output
            DAG._check_value(value, f'Output {name!r}')

    def _check_output_types(self):
        for name, output in self.outputs.items():
            if isinstance(output, OutputSpec) and output.type is None:
                raise ValueError(f'Output {name!r} has no explicit type')

    @staticmethod
    def _validate_mapped_ports(step):
        if step.map.get('group_by') is not None:
            return
        scatter = step.map.get('scatter')
        broadcast = step.map.get('broadcast')
        if scatter is None:
            raise ValueError(f'Mapped step {step.id!r} is missing map.scatter')
        if broadcast is None:
            raise ValueError(f'Mapped step {step.id!r} is missing map.broadcast')
        schema = step.input_schema or {}
        input_names = set(scatter) | set(broadcast)
        for name in schema:
            in_scatter = name in scatter
            in_broadcast = name in broadcast
            if in_scatter and in_broadcast:
                raise ValueError(
                    f'Mapped step {step.id!r}: input {name!r} appears in both '
                    f'scatter and broadcast'
                )
            if not in_scatter and not in_broadcast:
                raise ValueError(
                    f'Mapped step {step.id!r}: input {name!r} is neither '
                    f'scatter nor broadcast'
                )

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
            step = _step_from_dict(item, inputs, step_by_id)
            steps.append(step)
            step_by_id[step.id] = step
        for step, item in zip(steps, steps_data):
            step.bindings.update(_default_input_bindings(step, inputs, item))
            step.bindings.update({
                name: binding_from_dict(value, inputs, step_by_id)
                for name, value in item.get('bindings', {}).items()
            })
        outputs = {
            name: _output_spec_from_dict(value, inputs, step_by_id)
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


def _step_from_dict(item, inputs, step_by_id):
    step_outputs = item.get('outputs', {})
    step_run = item.get('run', {})
    step_inputs = item.get('inputs', {})
    map_data = item.get('map')
    return StepCall(
        id=item['id'],
        path=item['path'],
        bindings={},
        outputs=list(step_outputs.keys()),
        run={
            name: value
            for name, spec in step_run.items()
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
        **(_mapped_step_fields(map_data, item, inputs, step_by_id) if map_data is not None else {}),
    )



def _mapped_step_fields(map_data, item, inputs, step_by_id):
    source = None
    if map_data.get('source') is not None or map_data.get('fields') is not None:
        source = binding_from_dict(map_data['source'], inputs, step_by_id)
    return {
        'source': source,
        'map': map_data,
        'input_schema': item.get('input_schema'),
        'output_schema': item.get('output_schema'),
    }



def _default_input_bindings(step, inputs, item):
    return {
        name: inputs[name]
        for name in step.task.get('inputs', {}).keys()
        if name in inputs and name not in item.get('bindings', {})
    }



def _output_spec_from_dict(data, inputs, steps):
    if isinstance(data, dict) and 'value' in data:
        return OutputSpec(
            type=data.get('type'),
            desc=data.get('desc'),
            optional=data.get('optional', False),
            value=binding_from_dict(data['value'], inputs, steps),
        )
    return OutputSpec(type=None, value=binding_from_dict(data, inputs, steps))


def _run_param_to_dict(spec, value):
    data = dict(spec)
    default = data.pop('default', data.pop('value', None))
    if value is None:
        data['value'] = default
        return data
    if isinstance(value, Literal):
        data['value'] = value.value
        return data
    data['value'] = binding_to_dict(value)
    return data


def _run_value_from_dict(name, spec, defaults, inputs, steps):
    if not isinstance(spec, dict) or 'value' not in spec:
        return None
    value = spec['value']
    default = None
    if name in defaults and isinstance(defaults[name], dict):
        default = defaults[name].get('value')
    if isinstance(value, dict) and value.get('source') is not None:
        return binding_from_dict(value, inputs, steps)
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
