from dataclasses import dataclass, field
from typing import Dict, Optional



def _field_types(params):
    return {
        name: scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
        for name, param in params.items()
    }


@dataclass(frozen=True)
class WfType:
    pass


@dataclass(frozen=True)
class ScalarType(WfType):
    name: str


@dataclass(frozen=True)
class ArrayType(WfType):
    item: WfType


@dataclass(frozen=True)
class RecordType(WfType):
    fields: Dict[str, WfType] = field(default_factory=dict)
    open: bool = True


@dataclass(frozen=True)
class TableType(WfType):
    columns: Dict[str, WfType] = field(default_factory=dict)


@dataclass(frozen=True)
class FunctionType(WfType):
    input: WfType
    output: WfType


FILE = ScalarType('file')
STR = ScalarType('str')
INT = ScalarType('int')
FLOAT = ScalarType('float')
UNKNOWN = ScalarType('?')


_SCALARS = {
    'file': FILE,
    'str': STR,
    'int': INT,
    'float': FLOAT,
}


def scalar_from_name(name: Optional[str]) -> WfType:
    if name is None:
        return UNKNOWN
    if name.startswith('[') and name.endswith(']'):
        return ArrayType(scalar_from_name(name[1:-1]))
    return _SCALARS.get(name, ScalarType(name))


def from_task_signature(signature):
    return FunctionType(
        RecordType(_field_types(signature.inputs), open=True),
        RecordType(_field_types(signature.outputs), open=False),
    )


def is_record_type(value):
    return isinstance(value, RecordType)


def is_table_type(value):
    return isinstance(value, TableType)


def is_function_type(value):
    return isinstance(value, FunctionType)
