from dataclasses import dataclass, field
from typing import Dict, Optional


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
    inputs = {}
    for name, param in signature.inputs.items():
        typ = getattr(getattr(param, 'type', None), 'value', None)
        inputs[name] = scalar_from_name(typ)
    outputs = {}
    for name, param in signature.outputs.items():
        typ = getattr(getattr(param, 'type', None), 'value', None)
        outputs[name] = scalar_from_name(typ)
    return FunctionType(RecordType(inputs, open=True), RecordType(outputs, open=False))


def is_record_type(value):
    return isinstance(value, RecordType)


def is_table_type(value):
    return isinstance(value, TableType)


def is_function_type(value):
    return isinstance(value, FunctionType)
