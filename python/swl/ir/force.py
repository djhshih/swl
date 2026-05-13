from typing import Dict, List, Tuple

from swl.ir import node as ir
from swl.ir.dag import DAG, Field, ForcedFunction, Input, Literal, Merge, Record, TaskCall
from swl.ir.lower import Lowerer
from swl.semantic.task.type import signature_from_task
from swl.syntax.task import interpolation as interp
from swl.syntax.task.parser import Parser as TaskParser


_SENTINEL = object()


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
        self.variables = {}
        self.forced_variables = {}
        self.apply_cache = {}

    def force(self, node):
        self._assert_canonical(node)
        value = self.force_value(node, ForceEnv())
        value = self._force_root(value)
        self._refine_input_metadata()
        outputs = self._final_outputs(value)
        return DAG(dict(self.inputs), list(self.tasks), outputs)

    def _assert_canonical(self, node):
        if isinstance(node, ir.Chain):
            raise ValueError('Forcing expects normalized IR without Chain nodes')
        if isinstance(node, ir.Lambda):
            self._assert_canonical(node.body)
            return
        if isinstance(node, ir.Block):
            for bind in node.bindings:
                self._assert_canonical(bind)
            self._assert_canonical(node.result)
            return
        if isinstance(node, ir.Variable):
            self._assert_canonical(node.value)
            return
        if isinstance(node, ir.Apply):
            self._assert_canonical(node.function)
            self._assert_canonical(node.arg)
            return
        if isinstance(node, ir.Record):
            for value in node.fields.values():
                self._assert_canonical(value)
            return
        if isinstance(node, ir.Field):
            self._assert_canonical(node.record)
            return
        if isinstance(node, ir.Update):
            self._assert_canonical(node.left)
            self._assert_canonical(node.right)
            return

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

        if isinstance(node, ir.Ref):
            return self._force_ref(node, env)

        if isinstance(node, ir.Record):
            return Record({name: self.force_value(value, env) for name, value in node.fields.items()})

        if isinstance(node, ir.Field):
            source = self.force_value(node.record, env)
            projected = self._project_field(source, node.name)
            if projected is not _SENTINEL:
                return projected
            return Field(source, node.name)

        if isinstance(node, ir.Update):
            return self._merge_values(self.force_value(node.left, env), self.force_value(node.right, env))

        if isinstance(node, ir.Function):
            return ForcedFunction(node, None, node.signature)

        if isinstance(node, ir.Lambda):
            return ForcedFunction(node, None, node.signature)

        if isinstance(node, ir.Apply):
            return self._force_apply(node, env)

        if isinstance(node, ir.Block):
            local = ForceEnv(env)
            self._register_block(node.bindings)
            for bind in node.bindings:
                local.bind(bind.name, ir.Ref(bind.id, bind.name))
            return self.force_value(node.result, local)

        if isinstance(node, ir.Variable):
            self.variables[node.id] = node
            value = self.force_value(node.value, env)
            self.forced_variables[node.id] = value
            return value

        raise ValueError(f'Unsupported IR node during forcing: {type(node).__name__}')

    def _force_ref(self, ref, env):
        if ref.id in self.forced_variables:
            return self.forced_variables[ref.id]
        variable = self.variables[ref.id]
        value = self.force_value(variable.value, env)
        self.forced_variables[ref.id] = value
        return value

    def _force_apply(self, node, env):
        fn = self.force_value(node.function, env)
        arg = self.force_value(node.arg, env)
        key = self._apply_key(node.function, fn, arg)
        if key is not None and key in self.apply_cache:
            return self.apply_cache[key]
        result = self._apply(fn, arg)
        if key is not None:
            self.apply_cache[key] = result
        return result

    def _apply(self, fn, arg):
        if not isinstance(fn, ForcedFunction):
            raise ValueError(f'Cannot apply non-function value during forcing: {fn!r}')

        bound = self._merge_bound(fn.bound, arg)

        if isinstance(fn.function, ir.Function):
            if fn.function.kind == 'task':
                if not self._saturated_task(fn.function, bound):
                    return ForcedFunction(fn.function, bound, fn.signature)
                return self._emit_task_call(fn.function, bound)
            if fn.function.kind == 'workflow':
                if not self._saturated_signature(fn.function, bound):
                    return ForcedFunction(fn.function, bound, fn.signature)
                return self._force_workflow(fn.function, bound)

        if isinstance(fn.function, ir.Lambda):
            if not isinstance(bound, Record):
                return ForcedFunction(fn.function, bound, fn.signature)
            return self._force_lambda(fn.function, bound)

        return ForcedFunction(fn.function, bound, fn.signature)

    def _register_block(self, bindings):
        for bind in bindings:
            self.variables[bind.id] = bind

    def _forced_apply_key(self, value):
        if not isinstance(value, ForcedFunction):
            return None
        return _forced_function_key(value)

    def _apply_key(self, function_node, fn, arg):
        arg_key = _value_key(arg)
        forced_key = self._forced_apply_key(fn)
        if forced_key is not None:
            return ('apply-forced', forced_key, arg_key)
        if isinstance(function_node, ir.Ref):
            return ('apply-ref', function_node.id, arg_key)
        return None

    def _merge_bound(self, old, new):
        if old is None:
            return new
        return self._merge_values(old, new)

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

    def _project_field(self, source, name):
        if isinstance(source, Record):
            if name in source.fields:
                return source.fields[name]
            return _SENTINEL
        if isinstance(source, Merge):
            projected = self._project_field(source.right, name)
            if projected is not _SENTINEL:
                return projected
            return self._project_field(source.left, name)
        return _SENTINEL

    def _merge_values(self, left, right):
        if isinstance(left, Record) and isinstance(right, Record):
            fields = dict(left.fields)
            fields.update(right.fields)
            return Record(fields)
        return Merge(left, right)

    def _project_input(self, value, name):
        if isinstance(value, (Record, Merge)):
            projected = self._project_field(value, name)
            if projected is not _SENTINEL:
                return projected
            return Field(value, name)
        return value

    def _emit_task_call(self, function, bound):
        inputs = self._normalize_task_inputs(function, bound)
        key = self._task_call_key(function, inputs)
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
            return self._force_lambda(body, bound)
        return self.force_value(body, ForceEnv())

    def _task_call_key(self, function, inputs):
        return (
            function.path,
            tuple(sorted((name, _value_key(value)) for name, value in inputs.items())),
        )

    def _force_lambda(self, function, bound):
        local = ForceEnv()
        local.bind(function.param, bound)
        return self.force_value(function.body, local)

    def _task_definition(self, path):
        if path in self.task_defs:
            return self.task_defs[path]
        src = self.lowerer.checker._read_file(path)
        task = TaskParser().parse(src)
        signature = signature_from_task(task)
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
                    if section.kind.value == 'run':
                        semantic = signature.run[name]
                        target[name] = {
                            'type': semantic.type.value if semantic.type is not None else None,
                            'value': semantic.parsed_default,
                            'desc': param.desc,
                        }
                    else:
                        semantic = signature.inputs.get(name) if section.kind.value == 'in' else signature.outputs.get(name)
                        target[name] = {
                            'type': semantic.type.value if semantic and semantic.type is not None else param.type,
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

    def _refine_input_metadata(self):
        refined = {}
        for name, current in self.inputs.items():
            candidates = []
            for task in self.tasks:
                spec = task.task.get('inputs', {}).get(name)
                if spec is not None:
                    candidates.append(spec)
            best = current
            if candidates:
                typ = self._merge_input_type(name, current.type, candidates)
                desc = self._merge_input_desc(current.desc, candidates)
                best = Input(name, typ, desc)
            refined[name] = best
        self.inputs = refined

    def _merge_input_type(self, name, current, candidates):
        types = []
        if current is not None:
            types.append(current)
        types.extend(spec.get('type') for spec in candidates if spec.get('type') is not None)
        unique = list(dict.fromkeys(types))
        if not unique:
            return None
        if len(unique) > 1:
            raise ValueError(f'Conflicting input types during forcing: {name}: {unique}')
        return unique[0]

    def _merge_input_desc(self, current, candidates):
        descs = []
        if current:
            descs.append(current)
        descs.extend(spec.get('desc') for spec in candidates if spec.get('desc'))
        unique = list(dict.fromkeys(descs))
        if not unique:
            return None
        return ' / '.join(unique)

    def _task_dependencies(self, inputs):
        deps = set()
        for value in inputs.values():
            deps.update(self._value_dependencies(value))
        return deps

    def _walk_values(self, value):
        yield value
        if isinstance(value, Field):
            yield from self._walk_values(value.source)
        elif isinstance(value, Merge):
            yield from self._walk_values(value.left)
            yield from self._walk_values(value.right)
        elif isinstance(value, Record):
            for item in value.fields.values():
                yield from self._walk_values(item)

    def _value_dependencies(self, value):
        deps = set()
        for item in self._walk_values(value):
            if isinstance(item, Field) and isinstance(item.source, TaskCall):
                deps.add(item.source.id)
        return deps

    def _force_root(self, value):
        if not isinstance(value, ForcedFunction):
            return value
        signature = value.signature or self._function_signature(value.function)
        if signature is None:
            return value
        arg = Record({name: self._input(name, spec) for name, spec in signature.inputs.items()})
        return self._apply(value, arg)

    def _function_signature(self, function):
        if isinstance(function, ir.Function):
            return function.signature
        if isinstance(function, ir.Lambda):
            return function.signature
        return None

    def _input(self, name, spec=None):
        if name not in self.inputs:
            typ = None
            desc = None
            if spec is not None:
                typ = spec.type.value if getattr(spec, 'type', None) is not None else None
                desc = spec.desc
            self.inputs[name] = Input(name, typ, desc)
        return self.inputs[name]

    def _final_outputs(self, value):
        fields = self._collect_output_fields(value)
        if fields is not None:
            return fields
        return {'result': value}

    def _collect_output_fields(self, value):
        normalized = _normalize_output_value(value)
        if isinstance(normalized, Record):
            return dict(normalized.fields)
        return None


def _forced_function_key(value):
    function = value.function
    if isinstance(function, ir.Function):
        function_key = ('function', function.kind, function.path)
    elif isinstance(function, ir.Lambda):
        function_key = ('lambda', id(function))
    else:
        return None
    return (function_key, _value_key(value.bound))


def _flatten_value_terms(value):
    if isinstance(value, Merge):
        return _flatten_value_terms(value.left) + _flatten_value_terms(value.right)
    return [value]


def _record_fields_key(fields):
    return tuple(sorted((name, _value_key(value)) for name, value in fields.items()))


def _merge_key(value):
    terms = _flatten_value_terms(value)
    record_fields = {}
    others = []
    for term in terms:
        if isinstance(term, Record):
            record_fields.update(term.fields)
        else:
            others.append(_value_key(term))
    parts = []
    if record_fields:
        parts.append(('record', _record_fields_key(record_fields)))
    parts.extend(others)
    if len(parts) == 1:
        return parts[0]
    return ('merge', tuple(parts))


def _normalize_output_value(value):
    if isinstance(value, Merge):
        terms = _flatten_value_terms(value)
        record_fields = {}
        others = []
        for term in terms:
            normalized = _normalize_output_value(term)
            if isinstance(normalized, Record):
                record_fields.update(normalized.fields)
            else:
                others.append(normalized)
        if not others:
            return Record(record_fields)
        result = Record(record_fields) if record_fields else others[0]
        if record_fields:
            start = 0
        else:
            start = 1
        for term in others[start:]:
            result = Merge(result, term)
        return result
    if isinstance(value, Record):
        return Record({name: _normalize_output_value(item) for name, item in value.fields.items()})
    return value


def _value_key(value):
    if isinstance(value, Input):
        return ('input', value.name)
    if isinstance(value, Literal):
        return ('lit', value.value)
    if isinstance(value, Field):
        return ('field', _value_key(value.source), value.name)
    if isinstance(value, Merge):
        return _merge_key(value)
    if isinstance(value, Record):
        return ('record', _record_fields_key(value.fields))
    if isinstance(value, TaskCall):
        return ('task', value.path, value.id)
    return ('other', type(value).__name__, repr(value))


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


def force_file(path: str, files=None):
    forcer = Forcer(files=files)
    tree = forcer.lowerer.lower_file(path)
    return forcer.force(tree)
