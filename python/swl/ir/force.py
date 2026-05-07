import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from swl.ir import node as ir
from swl.ir.lower import Lowerer
from swl.syntax.task import interpolation as interp
from swl.syntax.task.parser import Parser as TaskParser


@dataclass(frozen=True)
class Input:
    name: str


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
    name: str
    path: str
    inputs: Dict[str, object]
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


@dataclass
class DAG:
    inputs: Dict[str, Input]
    tasks: List[TaskCall]
    outputs: Dict[str, object]

    def to_dict(self):
        return {
            'inputs': sorted(self.inputs.keys()),
            'tasks': [
                {
                    'id': task.id,
                    'name': task.name,
                    'path': task.path,
                    'inputs': {name: _value_to_dict(value) for name, value in task.inputs.items()},
                    'outputs': list(task.outputs),
                    'run': {name: _value_to_dict(value) for name, value in task.run.items()},
                    'task': task.task,
                    'deps': list(task.deps),
                }
                for task in self.tasks
            ],
            'outputs': {name: _value_to_dict(value) for name, value in self.outputs.items()},
        }

    @classmethod
    def from_dict(cls, data):
        inputs = {name: Input(name) for name in data.get('inputs', [])}
        tasks = []
        task_by_id = {}
        for item in data.get('tasks', []):
            task = TaskCall(
                id=item['id'],
                name=item['name'],
                path=item['path'],
                inputs={},
                outputs=list(item.get('outputs', [])),
                run={},
                task=item.get('task'),
                deps=list(item.get('deps', [])),
            )
            tasks.append(task)
            task_by_id[task.id] = task
        for task, item in zip(tasks, data.get('tasks', [])):
            task.inputs.update({
                name: _value_from_dict(value, inputs, task_by_id)
                for name, value in item.get('inputs', {}).items()
            })
            task.run.update({
                name: _value_from_dict(value, inputs, task_by_id)
                for name, value in item.get('run', {}).items()
            })
        outputs = {
            name: _value_from_dict(value, inputs, task_by_id)
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


class ForceEnv:
    def __init__(self, parent=None, values=None):
        self.parent = parent
        self.values = values or {}

    def bind(self, name, value):
        self.values[name] = value

    def lookup(self, name):
        if name in self.values:
            return self.values[name]
        if self.parent is not None:
            return self.parent.lookup(name)
        return None


class Forcer:
    def __init__(self, files=None):
        self.inputs = {}
        self.tasks = []
        self.task_cache = {}
        self.task_defs = {}
        self.call_counter = 0
        self.lowerer = Lowerer(files=files)

    def force(self, node):
        value = self.force_value(node, ForceEnv())
        value = self._force_root(value)
        outputs = self._final_outputs(value)
        return DAG(dict(self.inputs), list(self.tasks), outputs)

    def force_value(self, node, env):
        if isinstance(node, ir.Literal):
            return Literal(node.value)

        if isinstance(node, ir.Name):
            value = env.lookup(node.name)
            if value is not None:
                return value
            if node.name not in self.inputs:
                self.inputs[node.name] = Input(node.name)
            return self.inputs[node.name]

        if isinstance(node, ir.Record):
            return Record({name: self.force_value(value, env) for name, value in node.fields.items()})

        if isinstance(node, ir.Field):
            source = self.force_value(node.record, env)
            if isinstance(source, Record) and node.name in source.fields:
                return source.fields[node.name]
            return Field(source, node.name)

        if isinstance(node, ir.Update):
            return Merge(self.force_value(node.left, env), self.force_value(node.right, env))

        if isinstance(node, ir.Function):
            return ForcedFunction(node, None)

        if isinstance(node, ir.Lambda):
            return ForcedFunction(node, None)

        if isinstance(node, ir.Apply):
            fn = self.force_value(node.function, env)
            arg = self.force_value(node.arg, env)
            return self._apply(fn, arg)

        if isinstance(node, ir.Block):
            local = ForceEnv(env)
            for bind in node.bindings:
                local.bind(bind.name, self.force_value(bind.value, local))
            return self.force_value(node.result, local)

        if isinstance(node, ir.Chain):
            items = [self.force_value(item, env) for item in node.items]
            if not items:
                return Literal(None)
            result = items[0]
            for item in items[1:]:
                result = self._compose(result, item)
            return result

        if isinstance(node, ir.Bind):
            return self.force_value(node.value, env)

        return Literal(None)

    def _compose(self, left, right):
        return ForcedFunction(('chain', left, right), None)

    def _apply(self, fn, arg):
        if not isinstance(fn, ForcedFunction):
            return Literal(None)

        if isinstance(fn.function, tuple) and fn.function[0] == 'chain':
            left = self._apply(fn.function[1], arg)
            return self._apply(fn.function[2], left)

        bound = self._merge_bound(fn.bound, arg)

        if isinstance(fn.function, ir.Function):
            if fn.function.kind == 'task':
                if not self._saturated_task(fn.function, bound):
                    return ForcedFunction(fn.function, bound)
                return self._emit_task_call(fn.function, bound)
            if fn.function.kind == 'workflow':
                if not self._saturated_signature(fn.function, bound):
                    return ForcedFunction(fn.function, bound)
                return self._force_workflow(fn.function, bound)

        if isinstance(fn.function, ir.Lambda):
            if not isinstance(bound, Record):
                return ForcedFunction(fn.function, bound)
            local = ForceEnv()
            local.bind(fn.function.param, bound)
            return self.force_value(fn.function.body, local)

        return ForcedFunction(fn.function, bound)

    def _merge_bound(self, old, new):
        if old is None:
            return new
        return Merge(old, new)

    def _saturated_signature(self, function, bound):
        available = self._available_inputs(bound)
        return all(name in available for name in function.signature.inputs.keys())

    def _saturated_task(self, function, bound):
        available = self._available_inputs(bound)
        if len(available) == 1 and function.signature.inputs:
            first = next(iter(function.signature.inputs.keys()))
            if first in available:
                available = available
        return all(name in available for name in function.signature.inputs.keys())

    def _available_inputs(self, value):
        if isinstance(value, Input):
            return {value.name}
        if isinstance(value, Record):
            return set(value.fields.keys())
        if isinstance(value, Merge):
            return self._available_inputs(value.left).union(self._available_inputs(value.right))
        if isinstance(value, TaskCall):
            return set(value.outputs)
        if isinstance(value, Field):
            return set()
        if isinstance(value, Literal):
            return set()
        return set()

    def _normalize_task_inputs(self, function, bound):
        names = list(function.signature.inputs.keys())
        if not names:
            return {}
        if not self._looks_record_like(bound):
            return {names[0]: bound}
        return {name: self._project_input(bound, name) for name in names}

    def _looks_record_like(self, value):
        return isinstance(value, (Record, Merge))

    def _project_input(self, value, name):
        if isinstance(value, Record):
            if name in value.fields:
                return value.fields[name]
            return Field(value, name)
        if isinstance(value, Merge):
            right_fields = self._available_inputs(value.right)
            if name in right_fields:
                return self._project_input(value.right, name)
            left_fields = self._available_inputs(value.left)
            if name in left_fields:
                return self._project_input(value.left, name)
            return Field(value, name)
        return value

    def _emit_task_call(self, function, bound):
        inputs = self._normalize_task_inputs(function, bound)
        key = (function.name, tuple(sorted((name, _value_key(value)) for name, value in inputs.items())))
        if key in self.task_cache:
            return self.task_cache[key]
        self.call_counter += 1
        outputs = list(function.signature.outputs.keys())
        call = TaskCall(
            id=f't{self.call_counter}',
            name=function.name,
            path=function.path,
            inputs=inputs,
            outputs=outputs,
            run={},
            task=self._task_definition(function.path),
            deps=sorted(self._task_dependencies(inputs)),
        )
        self.tasks.append(call)
        result = Record({name: Field(call, name) for name in outputs})
        self.task_cache[key] = result
        return result

    def _force_workflow(self, function, bound):
        body = function.body
        if isinstance(body, ir.Lambda):
            local = ForceEnv()
            local.bind(body.param, bound)
            return self.force_value(body.body, local)
        return self.force_value(body, ForceEnv())

    def _task_definition(self, path):
        if path in self.task_defs:
            return self.task_defs[path]
        src = self.lowerer.checker._read_file(path)
        task = TaskParser().parse(src)
        inputs = {}
        outputs = {}
        run = {}
        for section in task.annotation.sections:
            target = None
            if section.kind.value == 'in':
                target = inputs
            elif section.kind.value == 'out':
                target = outputs
            elif section.kind.value == 'run':
                target = run
            for param in section.params:
                for name in param.names:
                    target[name] = {
                        'type': param.type,
                        'default': _interp_to_dict(param.default) if param.default is not None else None,
                        'desc': param.desc,
                    }
        definition = {
            'doc': task.annotation.doc,
            'body': task.body,
            'inputs': inputs,
            'outputs': outputs,
            'run': run,
        }
        self.task_defs[path] = definition
        return definition

    def _task_dependencies(self, inputs):
        deps = set()
        for value in inputs.values():
            deps.update(self._value_dependencies(value))
        return deps

    def _value_dependencies(self, value):
        if isinstance(value, Field) and isinstance(value.source, TaskCall):
            return {value.source.id}
        if isinstance(value, Field):
            return self._value_dependencies(value.source)
        if isinstance(value, Merge):
            return self._value_dependencies(value.left).union(self._value_dependencies(value.right))
        if isinstance(value, Record):
            deps = set()
            for item in value.fields.values():
                deps.update(self._value_dependencies(item))
            return deps
        return set()

    def _force_root(self, value):
        if not isinstance(value, ForcedFunction):
            return value
        signature = self._function_signature(value.function)
        if signature is None:
            return value
        arg = Record({name: self._input(name) for name in signature.inputs.keys()})
        return self._apply(value, arg)

    def _function_signature(self, function):
        if isinstance(function, ir.Function):
            return function.signature
        if isinstance(function, ir.Lambda):
            return function.signature
        return None

    def _input(self, name):
        if name not in self.inputs:
            self.inputs[name] = Input(name)
        return self.inputs[name]

    def _final_outputs(self, value):
        if isinstance(value, Record):
            return dict(value.fields)
        return {'result': value}


def _value_key(value):
    if isinstance(value, Input):
        return ('input', value.name)
    if isinstance(value, Literal):
        return ('lit', value.value)
    if isinstance(value, Field):
        return ('field', _value_key(value.source), value.name)
    if isinstance(value, Merge):
        return ('merge', _value_key(value.left), _value_key(value.right))
    if isinstance(value, Record):
        return ('record', tuple(sorted((name, _value_key(v)) for name, v in value.fields.items())))
    if isinstance(value, TaskCall):
        return ('task', value.id)
    return ('other', repr(value))


def _value_to_dict(value):
    if isinstance(value, Input):
        return {'kind': 'input', 'name': value.name}
    if isinstance(value, Literal):
        return {'kind': 'literal', 'value': value.value}
    if isinstance(value, Field):
        if isinstance(value.source, TaskCall):
            return {'kind': 'task_output', 'task': value.source.id, 'name': value.name}
        return {'kind': 'field', 'source': _value_to_dict(value.source), 'name': value.name}
    if isinstance(value, Merge):
        return {'kind': 'merge', 'left': _value_to_dict(value.left), 'right': _value_to_dict(value.right)}
    if isinstance(value, Record):
        return {'kind': 'record', 'fields': {name: _value_to_dict(v) for name, v in value.fields.items()}}
    if isinstance(value, TaskCall):
        return {'kind': 'task', 'id': value.id}
    if isinstance(value, ForcedFunction):
        return {'kind': 'function'}
    return {'kind': 'unknown', 'repr': repr(value)}


def _interp_to_dict(value):
    if isinstance(value, interp.Word):
        return {'kind': 'word', 'parts': [_interp_to_dict(part) for part in value.parts]}
    if isinstance(value, interp.Literal):
        return {'kind': 'literal', 'text': value.text}
    if isinstance(value, interp.Var):
        return {'kind': 'var', 'name': value.name}
    if isinstance(value, interp.Expr):
        return {'kind': 'expr', 'text': value.text}
    return None


def _value_from_dict(data, inputs, tasks):
    kind = data.get('kind')
    if kind == 'input':
        return inputs[data['name']]
    if kind == 'literal':
        return Literal(data.get('value'))
    if kind == 'task_output':
        return Field(tasks[data['task']], data['name'])
    if kind == 'field':
        return Field(_value_from_dict(data['source'], inputs, tasks), data['name'])
    if kind == 'merge':
        return Merge(
            _value_from_dict(data['left'], inputs, tasks),
            _value_from_dict(data['right'], inputs, tasks),
        )
    if kind == 'record':
        return Record({
            name: _value_from_dict(value, inputs, tasks)
            for name, value in data.get('fields', {}).items()
        })
    if kind == 'task':
        return tasks[data['id']]
    if kind == 'function':
        return {'kind': 'function'}
    return {'kind': 'unknown', 'repr': data.get('repr')}


def force_file(path: str, files=None):
    forcer = Forcer(files=files)
    tree = forcer.lowerer.lower_file(path)
    return forcer.force(tree)
