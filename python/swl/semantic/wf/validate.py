from swl.semantic.wf import type as wf_type


class WorkflowInputValidationError(ValueError):
    pass


class WorkflowInputValidator:
    def validate(self, workflow_check, inputs):
        workflow_type = getattr(workflow_check, 'workflow_type', None)
        if not isinstance(workflow_type, wf_type.FunctionType):
            raise WorkflowInputValidationError('Workflow must have a function type before inputs can be validated')
        return self._validate_value(workflow_type.input, inputs, path='input')

    def _validate_value(self, expected, value, path):
        if isinstance(expected, wf_type.RecordType):
            return self._validate_record(expected, value, path)
        if isinstance(expected, wf_type.TableType):
            return self._validate_table(expected, value, path)
        if isinstance(expected, wf_type.ArrayType):
            return self._validate_array(expected, value, path)
        if isinstance(expected, wf_type.ScalarType):
            return self._validate_scalar(expected, value, path)
        return value

    def _validate_record(self, expected, value, path):
        if not isinstance(value, dict):
            raise WorkflowInputValidationError(f'{path} must be a record/object')
        for name, field_type in expected.fields.items():
            if name not in value:
                raise WorkflowInputValidationError(f'{path}.{name} is required')
            self._validate_value(field_type, value[name], f'{path}.{name}')
        return value

    def _validate_table(self, expected, value, path):
        if not isinstance(value, dict):
            raise WorkflowInputValidationError(f'{path} must be a table/object')
        # Current checker may expose a placeholder single-name root table parameter
        # (for example TableType(columns={'xs': ?})) when concrete columns are not
        # yet propagated. Accept either the wrapped form {xs: {...}} or a direct
        # concrete column object for forward compatibility.
        if len(expected.columns) == 1:
            only_name, only_type = next(iter(expected.columns.items()))
            if only_name in value and isinstance(value[only_name], dict):
                return self._validate_unknown_table(value[only_name], f'{path}.{only_name}')
            if only_type == wf_type.UNKNOWN:
                return self._validate_unknown_table(value, path)
        lengths = set()
        for name, column_type in expected.columns.items():
            if name not in value:
                raise WorkflowInputValidationError(f'{path}.{name} is required')
            column = value[name]
            if not isinstance(column, list):
                raise WorkflowInputValidationError(f'{path}.{name} must be an array column')
            lengths.add(len(column))
            self._validate_array(wf_type.ArrayType(column_type), column, f'{path}.{name}')
        if len(lengths) > 1:
            raise WorkflowInputValidationError(f'{path} table columns must have equal length')
        return value

    def _validate_unknown_table(self, value, path):
        if not isinstance(value, dict):
            raise WorkflowInputValidationError(f'{path} must be a table/object')
        lengths = set()
        for name, column in value.items():
            if not isinstance(column, list):
                raise WorkflowInputValidationError(f'{path}.{name} must be an array column')
            lengths.add(len(column))
        if len(lengths) > 1:
            raise WorkflowInputValidationError(f'{path} table columns must have equal length')
        return value

    def _validate_array(self, expected, value, path):
        if not isinstance(value, list):
            raise WorkflowInputValidationError(f'{path} must be an array')
        for index, item in enumerate(value):
            self._validate_value(expected.item, item, f'{path}[{index}]')
        return value

    def _validate_scalar(self, expected, value, path):
        name = expected.name
        if name == '?':
            return value
        if name == 'file' and isinstance(value, str):
            return value
        if name == 'str' and isinstance(value, str):
            return value
        if name == 'int' and isinstance(value, int) and not isinstance(value, bool):
            return value
        if name == 'float' and isinstance(value, (int, float)) and not isinstance(value, bool):
            return value
        raise WorkflowInputValidationError(f'{path} must have type {name}')


def validate_workflow_inputs(workflow_check, inputs):
    return WorkflowInputValidator().validate(workflow_check, inputs)
