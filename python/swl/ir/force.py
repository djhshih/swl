from typing import Dict, List, Tuple

from swl.ir import node as ir
from swl.ir.dag import DAG, ArrayField, Field, ForcedFunction, Input, Literal, Merge, Record, StepCall, MappedStep, MappedValue
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
        self.steps = []
        self.step_cache = {}
        self.task_defs = {}
        self.call_counter = 0
        self.task_name_counts = {}
        self.lowerer = Lowerer(files=files)
        self.variables = {}
        self.forced_variables = {}
        self.apply_cache = {}

    def force(self, node):
        self._assert_canonical(node)
        value = self.force_value(node, ForceEnv())
        value = self._force_root(value)
        return self._finalize_dag(value)

    def _finalize_dag(self, value):
        self._refine_input_metadata()
        outputs = self._final_outputs(value)
        self._prune_unused_inputs(value, outputs)
        return DAG(dict(self.inputs), list(self.steps), outputs)

    def _assert_canonical(self, node):
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
        if isinstance(node, ir.ArrayField):
            self._assert_canonical(node.record_array)
            return
        if isinstance(node, ir.Map):
            self._assert_canonical(node.function)
            self._assert_canonical(node.arg)
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
            if isinstance(source, MappedStep):
                if node.name not in source.outputs:
                    raise ValueError(f'Missing field on mapped result: {node.name}')
                return ArrayField(source, node.name)
            if isinstance(source, MappedValue):
                projected = self._project_field(source.element, node.name)
                if projected is _SENTINEL and isinstance(source.element, Field) and source.element.name == node.name:
                    projected = source.element
                if projected is _SENTINEL and isinstance(source.element, Record) and node.name in source.element.fields:
                    projected = source.element.fields[node.name]
                if projected is _SENTINEL:
                    raise ValueError(f'Missing field on mapped result: {node.name}')
                return MappedValue(source.source, projected)
            return Field(source, node.name)

        if isinstance(node, ir.ArrayField):
            source = self.force_value(node.record_array, env)
            if isinstance(source, MappedStep):
                if node.name not in source.outputs:
                    raise ValueError(f'Missing field on mapped result: {node.name}')
                return ArrayField(source, node.name)
            return ArrayField(source, node.name)

        if isinstance(node, ir.Update):
            return self._merge_values(self.force_value(node.left, env), self.force_value(node.right, env))

        if isinstance(node, ir.Function):
            return ForcedFunction(node, None, node.signature)

        if isinstance(node, ir.Lambda):
            return ForcedFunction(node, None, node.signature)

        if isinstance(node, ir.Apply):
            return self._force_apply(node, env)

        if isinstance(node, ir.Map):
            fn = self.force_value(node.function, env)
            arg = self.force_value(node.arg, env)
            return self._force_map(fn, arg)

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
        if isinstance(arg, (MappedStep, MappedValue)):
            return self._apply_mapped(fn, arg)
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
                return self._emit_workflow_call(fn.function, bound)

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
        if self._is_opaque_record_carrier(bound):
            return True
        available = self._available_inputs(bound)
        return all(name in available for name in function.signature.inputs.keys())

    def _saturated_task(self, function, bound):
        if self._is_opaque_record_carrier(bound):
            return True
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
        if isinstance(value, StepCall):
            return set(value.outputs)
        if isinstance(value, MappedStep):
            return set(value.outputs)
        if isinstance(value, MappedValue):
            projected = self._collect_output_fields(value.element)
            return set(projected.keys()) if projected is not None else set()
        if isinstance(value, Field):
            return set()
        if isinstance(value, ArrayField):
            return set()
        if isinstance(value, Literal):
            return set()
        return set()

    # task argument shaping -------------------------------------------------

    def _normalize_task_inputs(self, function, bound):
        names = list(function.signature.inputs.keys())
        if not names:
            return {}
        if not self._can_project_record_fields(bound):
            return {names[0]: bound}
        return {name: self._project_input(bound, name) for name in names}

    def _looks_record_like(self, value):
        return isinstance(value, (Record, Merge))

    def _can_project_record_fields(self, value):
        return isinstance(value, (Input, Field, Record, Merge))

    def _is_opaque_record_carrier(self, value):
        if isinstance(value, (Input, Field)):
            return True
        if isinstance(value, Record) and len(value.fields) == 1:
            only = next(iter(value.fields.values()))
            return isinstance(only, (Input, Field))
        return False

    def _project_field(self, source, name):
        if isinstance(source, Input):
            return Field(source, name)
        if isinstance(source, Field):
            return Field(source, name)
        if isinstance(source, Record):
            if name in source.fields:
                return source.fields[name]
            if len(source.fields) == 1:
                only = next(iter(source.fields.values()))
                if isinstance(only, (Input, Field)):
                    return Field(only, name)
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
        if self._can_project_record_fields(value):
            projected = self._project_field(value, name)
            if projected is not _SENTINEL:
                return projected
            return Field(value, name)
        return value

    def _emit_task_call(self, function, bound):
        inputs = self._normalize_task_inputs(function, bound)
        key = self._step_call_key(function, inputs)
        if key in self.step_cache:
            return self.step_cache[key]
        self.call_counter += 1
        outputs = list(function.signature.outputs.keys())
        call_id = self._task_id(function.name)
        call = StepCall(
            id=call_id,
            path=function.path,
            bindings=inputs,
            outputs=outputs,
            run=self._normalize_task_run(function, bound),
            task=self._task_definition(function.path),
            deps=sorted(self._step_dependencies(inputs)),
        )
        self.steps.append(call)
        result = Record({name: Field(call, name) for name in outputs})
        self.step_cache[key] = result
        return result

    def _force_map(self, fn, arg):
        if not isinstance(fn, ForcedFunction):
            raise ValueError(f'Cannot map non-function value during forcing: {fn!r}')
        source = arg
        if isinstance(arg, Record) and 'xs' in arg.fields:
            source = arg.fields['xs']
        return self._apply_mapped(fn, source)

    def _apply_mapped(self, fn, source):
        if not isinstance(fn, ForcedFunction):
            raise ValueError(f'Cannot map non-function value during forcing: {fn!r}')
        if self._is_batch_function(fn):
            raise ValueError('map on batch workflow is not supported during forcing')
        if isinstance(fn.function, ir.Function) and fn.function.kind in ('task', 'workflow'):
            return self._emit_mapped_step(fn, source)
        raise ValueError(f'map requires normalized executable callable during forcing: {fn.function!r}')

    def _emit_mapped_step(self, fn, source):
        target = fn.function
        outputs = list(target.signature.outputs.keys())
        step_id = self._task_id(target.name)
        step = MappedStep(
            id=step_id,
            path=target.path,
            source=source,
            outputs=outputs,
            task=self._task_definition(target.path) if target.kind == 'task' else self._workflow_definition(target),
            deps=sorted(self._step_dependencies({'source': source})),
            type=target.kind,
            map={'source': _binding_to_public_dict(source)},
        )
        self.steps.append(step)
        return step

    def _task_id(self, name):
        count = self.task_name_counts.get(name, 0) + 1
        self.task_name_counts[name] = count
        if count == 1:
            return name
        return f'{name}_{count}'

    def _emit_workflow_call(self, function, bound):
        inputs = self._normalize_task_inputs(function, bound)
        key = self._step_call_key(function, inputs)
        if key in self.step_cache:
            return self.step_cache[key]
        outputs = list(function.signature.outputs.keys())
        call_id = self._task_id(function.name)
        call = StepCall(
            id=call_id,
            path=function.path,
            bindings=inputs,
            outputs=outputs,
            run={},
            task=self._workflow_definition(function),
            deps=sorted(self._step_dependencies(inputs)),
            type='workflow',
        )
        self.steps.append(call)
        result = Record({name: Field(call, name) for name in outputs})
        self.step_cache[key] = result
        return result

    def _force_workflow(self, function, bound):
        body = function.body
        if isinstance(body, ir.Lambda):
            return self._force_lambda(body, bound)
        return self.force_value(body, ForceEnv())

    def _mapped_element_input(self, source):
        if isinstance(source, Input):
            return Record({source.name: Input(source.name)})
        if isinstance(source, ArrayField):
            return Field(source.source, source.name)
        if isinstance(source, MappedValue):
            return source.element
        if isinstance(source, MappedStep):
            return Record({name: Field(source, name) for name in source.outputs})
        return source

    def _is_batch_function(self, fn):
        signature = self._forced_signature(fn)
        if signature is None:
            return False
        input_names = list(signature.inputs.keys())
        return input_names == ['xs']

    def _step_call_key(self, function, inputs):
        return (
            function.path,
            tuple(sorted((name, _value_key(value)) for name, value in inputs.items())),
        )

    def _force_lambda(self, function, bound):
        local = ForceEnv()
        local.bind(function.param, bound)
        return self.force_value(function.body, local)

    def _normalize_task_run(self, function, bound):
        names = list(function.signature.run.keys())
        if not names or not self._can_project_record_fields(bound):
            return {}
        return {
            name: self._project_input(bound, name)
            for name in names
            if self._project_field(bound, name) is not _SENTINEL
        }

    def _task_definition(self, path):
        if path in self.task_defs:
            return self.task_defs[path]
        src = self.lowerer.checker._read_file(path)
        task = TaskParser().parse(src)
        signature = signature_from_task(task)
        definition = {
            'doc': task.annotation.doc,
            'body': task.body,
            'inputs': self._build_task_section(task, signature, 'in'),
            'outputs': self._build_task_section(task, signature, 'out'),
            'run': self._build_task_section(task, signature, 'run'),
        }
        self.task_defs[path] = definition
        return definition

    def _materialize_workflow_dag(self, function):
        forcer = Forcer(files=self.lowerer.checker.files)
        signature = function.signature
        if signature is None:
            return forcer.force(function.body)
        arg = Record({name: forcer._input(name, spec) for name, spec in signature.inputs.items()})
        if isinstance(function.body, ir.Lambda):
            value = forcer._apply(ForcedFunction(function.body, None, signature), arg)
        else:
            value = forcer.force_value(function.body, ForceEnv())
            if isinstance(value, ForcedFunction):
                value = forcer._apply(value, arg)
        return forcer._finalize_dag(value)

    def _workflow_definition(self, function):
        key = ('workflow', function.path)
        if key in self.task_defs:
            return self.task_defs[key]
        body_dag = self._materialize_workflow_dag(function).to_dict()
        definition = {
            'class': 'Workflow',
            'dag': body_dag,
            'inputs': {name: {'type': self._normalize_swl_type(spec.type.value if getattr(spec, 'type', None) is not None else None), 'desc': spec.desc} for name, spec in function.signature.inputs.items()},
            'outputs': {name: {'type': self._normalize_swl_type(spec.type.value if getattr(spec, 'type', None) is not None else None)} for name, spec in function.signature.outputs.items()},
            'body': '',
            'run': {},
        }
        self.task_defs[key] = definition
        return definition

    def _build_task_section(self, task, signature, kind):
        built = {}
        for section in task.annotation.sections:
            if section.kind.value != kind:
                continue
            for param in section.params:
                for name in param.names:
                    built[name] = self._build_task_param(signature, kind, name, param)
        return built

    def _normalize_swl_type(self, value):
        if value == 'File':
            return 'file'
        if value == 'string':
            return 'str'
        return value

    def _build_task_param(self, signature, kind, name, param):
        if kind == 'run':
            semantic = signature.run[name]
            built = {
                'type': semantic.type.value if semantic.type is not None else None,
                'value': semantic.parsed_default,
            }
            if param.desc is not None:
                built['desc'] = param.desc
            return built
        semantic = signature.inputs.get(name) if kind == 'in' else signature.outputs.get(name)
        built = {
            'type': semantic.type.value if semantic and semantic.type is not None else param.type,
        }
        if param.desc is not None:
            built['desc'] = param.desc
        if param.default is not None:
            built['default'] = _interp_to_dict(param.default)
        return built

    # dag finalization ------------------------------------------------------

    def _refine_input_metadata(self):
        refined = {}
        for name, current in self.inputs.items():
            candidates = []
            for step in self.steps:
                spec = step.task.get('inputs', {}).get(name)
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

    def _prune_unused_inputs(self, root_value, outputs):
        used = set()
        if isinstance(root_value, ForcedFunction):
            signature = self._forced_signature(root_value)
            if signature is not None:
                used.update(signature.inputs.keys())
        for step in self.steps:
            for value in step.bindings.values():
                for item in self._walk_values(value):
                    if isinstance(item, Input):
                        used.add(item.name)
            mapped_source = getattr(step, 'source', None)
            if mapped_source is not None:
                for item in self._walk_values(mapped_source):
                    if isinstance(item, Input):
                        used.add(item.name)
        for value in outputs.values():
            for item in self._walk_values(value):
                if isinstance(item, Input):
                    used.add(item.name)
        for step in self.steps:
            mapped = getattr(step, 'map', None)
            if mapped is not None:
                source = mapped.get('source', {})
                if source.get('source') == 'input' and source.get('name') not in self.inputs:
                    self.inputs[source['name']] = Input(source['name'])
                if source.get('source') == 'input':
                    used.add(source['name'])
        self.inputs = {name: value for name, value in self.inputs.items() if name in used}

    # dependency and reachability helpers ----------------------------------

    def _step_dependencies(self, inputs):
        deps = set()
        for value in inputs.values():
            deps.update(self._value_dependencies(value))
        return deps

    def _walk_values(self, value):
        yield value
        if isinstance(value, Field):
            yield from self._walk_values(value.source)
        elif isinstance(value, ArrayField):
            yield from self._walk_values(value.source)
        elif isinstance(value, MappedValue):
            yield from self._walk_values(value.source)
            yield from self._walk_values(value.element)
        elif isinstance(value, Merge):
            yield from self._walk_values(value.left)
            yield from self._walk_values(value.right)
        elif isinstance(value, Record):
            for item in value.fields.values():
                yield from self._walk_values(item)

    def _value_dependencies(self, value):
        deps = set()
        for item in self._walk_values(value):
            if isinstance(item, Field) and isinstance(item.source, (StepCall, MappedStep)):
                deps.add(item.source.id)
            if isinstance(item, ArrayField) and isinstance(item.source, (StepCall, MappedStep)):
                deps.add(item.source.id)
        return deps

    # root forcing and signature helpers -----------------------------------

    def _force_root(self, value):
        if not isinstance(value, ForcedFunction):
            return value
        signature = self._forced_signature(value)
        if signature is None:
            return value
        arg = Record({name: self._input(name, spec) for name, spec in signature.inputs.items()})
        return self._apply(value, arg)

    def _forced_signature(self, value):
        signature = value.signature or self._function_signature(value.function)
        if signature is None:
            return None
        if value.bound is None:
            return signature
        available = self._available_inputs(value.bound)
        inputs = {
            name: param
            for name, param in signature.inputs.items()
            if name not in available
        }
        return type(signature)(inputs, dict(signature.outputs), dict(signature.run))

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
    if isinstance(value, Field):
        source = _normalize_output_value(value.source)
        if isinstance(source, Record) and len(source.fields) == 1:
            only = next(iter(source.fields.values()))
            if isinstance(only, Input):
                return Field(only, value.name)
        return Field(source, value.name)
    return value


def _value_key(value):
    if isinstance(value, Input):
        return ('input', value.name)
    if isinstance(value, Literal):
        return ('lit', value.value)
    if isinstance(value, Field):
        return ('field', _value_key(value.source), value.name)
    if isinstance(value, ArrayField):
        return ('array_field', _value_key(value.source), value.name)
    if isinstance(value, Merge):
        return _merge_key(value)
    if isinstance(value, Record):
        return ('record', _record_fields_key(value.fields))
    if isinstance(value, (StepCall, MappedStep)):
        return ('step', value.path, value.id)
    if isinstance(value, MappedValue):
        return ('mapped', _value_key(value.source), _value_key(value.element))
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


def _binding_to_public_dict(value):
    if isinstance(value, Input):
        return {'source': 'input', 'name': value.name}
    if isinstance(value, Field) and isinstance(value.source, (StepCall, MappedStep)):
        return {'source': 'step', 'step': value.source.id, 'output': value.name}
    if isinstance(value, ArrayField) and isinstance(value.source, (StepCall, MappedStep)):
        return {'source': 'array_field', 'step': value.source.id, 'output': value.name}
    return {'source': 'value'}


def force_file(path: str, files=None):
    forcer = Forcer(files=files)
    tree = forcer.lowerer.lower_file(path)
    return forcer.force(tree)
