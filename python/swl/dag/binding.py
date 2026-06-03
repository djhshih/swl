def is_instance(value, class_name):
    return value.__class__.__name__ == class_name


def binding_to_dict(value):
    cls = value.__class__.__name__
    if cls == 'Input':
        return {'source': 'input', 'name': value.name}
    if cls == 'Literal':
        return {'source': 'literal', 'value': value.value}
    if cls == 'Field':
        if is_instance(value.source, 'StepCall'):
            return {'source': 'step_output', 'step': value.source.id, 'output': value.name}
        return {'source': 'field', 'field': value.name, 'value': binding_to_dict(value.source)}
    if cls == 'Merge':
        return {'source': 'merge', 'left': binding_to_dict(value.left), 'right': binding_to_dict(value.right)}
    if cls == 'Record':
        return {'source': 'record', 'fields': {name: binding_to_dict(v) for name, v in value.fields.items()}}
    if cls == 'TableSource':
        return {'source': 'table', 'name': value.name, 'columns': {name: binding_to_dict(v) for name, v in value.columns.items()}}
    if cls == 'StepCall':
        return {'source': 'step_call', 'step': value.id}
    if cls == 'ForcedFunction':
        return {'source': 'function'}
    raise ValueError(f'Unsupported value for binding serialization: {cls}')


def binding_from_dict(data, inputs, steps):
    if not isinstance(data, dict):
        from swl.dag.node import Literal
        return Literal(data)
    source = data.get('source')
    if source == 'input':
        from swl.dag.node import Input
        name = data['name']
        return inputs.get(name, Input(name, optional=False))
    if source == 'value':
        return binding_from_dict(data['value'], inputs, steps)
    if source == 'literal':
        from swl.dag.node import Literal
        return Literal(data.get('value'))
    if source == 'field':
        from swl.dag.node import Field
        return Field(binding_from_dict(data['value'], inputs, steps), data['field'])
    if source == 'merge':
        from swl.dag.node import Merge
        return Merge(
            binding_from_dict(data['left'], inputs, steps),
            binding_from_dict(data['right'], inputs, steps),
        )
    if source == 'record':
        from swl.dag.node import Record
        return Record({
            name: binding_from_dict(value, inputs, steps)
            for name, value in data.get('fields', {}).items()
        })
    if source == 'table':
        from swl.dag.node import TableSource
        return TableSource(
            data['name'],
            {name: binding_from_dict(value, inputs, steps) for name, value in data.get('columns', {}).items()},
        )
    if source == 'step_output':
        from swl.dag.node import Field
        step = steps[data['step']]
        return Field(step, data['output'])
    if source == 'step_call':
        return steps[data['step']]
    if source == 'function':
        return {'kind': 'function'}
    if 'step' in data and 'output' in data:
        from swl.dag.node import Field
        step = steps[data['step']]
        return Field(step, data['output'])
    raise ValueError(f'Unsupported binding source during deserialization: {source!r}')
