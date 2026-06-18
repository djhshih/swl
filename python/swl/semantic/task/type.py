import math
import re
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Set

from swl.syntax.task import interpolation as interp
from swl.syntax.task import node as task_node


class TypeKind(Enum):
    FILE = 'file'
    STR = 'str'
    INT = 'int'
    FLOAT = 'float'
    MEMORY = 'memory'
    TIME = 'time'
    FILE_OPT = 'file?'
    STR_OPT = 'str?'
    INT_OPT = 'int?'
    FLOAT_OPT = 'float?'
    ARRAY_FILE = '[file]'
    ARRAY_STR = '[str]'
    ARRAY_INT = '[int]'
    ARRAY_FLOAT = '[float]'


TYPE_MAP = {kind.value: kind for kind in TypeKind}


OPTIONAL_TYPES = {
    TypeKind.FILE_OPT,
    TypeKind.STR_OPT,
    TypeKind.INT_OPT,
    TypeKind.FLOAT_OPT,
}


ARRAY_TYPES = {
    TypeKind.ARRAY_FILE,
    TypeKind.ARRAY_STR,
    TypeKind.ARRAY_INT,
    TypeKind.ARRAY_FLOAT,
}


BASE_TYPE_MAP = {
    TypeKind.FILE_OPT: TypeKind.FILE,
    TypeKind.STR_OPT: TypeKind.STR,
    TypeKind.INT_OPT: TypeKind.INT,
    TypeKind.FLOAT_OPT: TypeKind.FLOAT,
    TypeKind.ARRAY_FILE: TypeKind.FILE,
    TypeKind.ARRAY_STR: TypeKind.STR,
    TypeKind.ARRAY_INT: TypeKind.INT,
    TypeKind.ARRAY_FLOAT: TypeKind.FLOAT,
}


def parse_type(type_str: str) -> TypeKind:
    '''Parse a type string into a TypeKind.'''
    type_str = type_str.strip()
    if type_str not in TYPE_MAP:
        raise ValueError(f'Unknown type: {type_str}')
    return TYPE_MAP[type_str]


def is_optional(typ: TypeKind) -> bool:
    return typ in OPTIONAL_TYPES


def is_array(typ: TypeKind) -> bool:
    return typ in ARRAY_TYPES


def base_type(typ: TypeKind) -> TypeKind:
    return BASE_TYPE_MAP.get(typ, typ)


COMPATIBLE_PAIRS = {
    (TypeKind.FILE, TypeKind.FILE): True,
    (TypeKind.FILE, TypeKind.FILE_OPT): True,
    (TypeKind.STR, TypeKind.STR): True,
    (TypeKind.STR, TypeKind.STR_OPT): True,
    (TypeKind.INT, TypeKind.INT): True,
    (TypeKind.INT, TypeKind.INT_OPT): True,
    (TypeKind.FLOAT, TypeKind.FLOAT): True,
    (TypeKind.FLOAT, TypeKind.FLOAT_OPT): True,
    (TypeKind.MEMORY, TypeKind.MEMORY): True,
    (TypeKind.TIME, TypeKind.TIME): True,
}


def types_compatible(output_type: TypeKind, input_type: TypeKind) -> bool:
    if (output_type, input_type) in COMPATIBLE_PAIRS:
        return True

    output_base = base_type(output_type)
    input_base = base_type(input_type)

    if is_optional(input_type) and output_base == input_base:
        return True

    return False


@dataclass
class Param:
    '''Semantic task parameter.'''

    name: str
    type: TypeKind | None = None
    default: str | None = None
    desc: str | None = None
    parsed_default: object = None


@dataclass
class TaskSignature:
    '''Type signature for a task.'''

    inputs: Dict[str, Param]
    outputs: Dict[str, Param]
    run: Dict[str, Param] | None = None

    def __post_init__(self):
        if self.run is None:
            self.run = {}

    def input_names(self) -> Set[str]:
        return set(self.inputs.keys())

    def output_names(self) -> Set[str]:
        return set(self.outputs.keys())


class TypeChecker:
    '''Type checker for workflows.'''

    def __init__(self):
        self.signatures: Dict[str, TaskSignature] = {}

    def add_task(self, name: str, signature: TaskSignature):
        self.signatures[name] = signature

    def check_chain(self, left: str, right: str) -> List[str]:
        errors = []

        if left not in self.signatures:
            errors.append(f'Unknown task: {left}')
            return errors

        if right not in self.signatures:
            errors.append(f'Unknown task: {right}')
            return errors

        left_sig = self.signatures[left]
        right_sig = self.signatures[right]

        for out_name, out_param in left_sig.outputs.items():
            if out_name in right_sig.inputs:
                in_param = right_sig.inputs[out_name]
                if not types_compatible(out_param.type, in_param.type):
                    errors.append(
                        f'Type mismatch for "{out_name}": '
                        f'{out_param.type.value} -> {in_param.type.value}'
                    )

        return errors

def signature_from_task(task: task_node.Task) -> TaskSignature:
    '''Build semantic task signature from parsed task syntax.'''
    inputs = {}
    outputs = {}
    run = {}
    in_names = set()
    out_names = set()
    run_names = set()

    for section in task.annotation.sections:
        for parsed in section.params:
            typ = parse_type(parsed.type) if parsed.type else None
            for name in parsed.names:
                if section.kind == task_node.SectionType.IN:
                    if name in in_names:
                        raise ValueError(f'Duplicate input parameter: {name}')
                    if name in run_names:
                        raise ValueError(f'Duplicate parameter: {name}')
                    in_names.add(name)
                elif section.kind == task_node.SectionType.OUT:
                    if name in out_names:
                        raise ValueError(f'Duplicate output parameter: {name}')
                    if name in run_names:
                        raise ValueError(f'Duplicate parameter: {name}')
                    out_names.add(name)
                elif section.kind == task_node.SectionType.RUN:
                    if name in run_names:
                        raise ValueError(f'Duplicate run parameter: {name}')
                    if name in in_names or name in out_names:
                        raise ValueError(f'Duplicate parameter: {name}')
                    run_names.add(name)
                else:
                    raise ValueError(f'Unrecognized section kind: {section.kind}')
                parsed_default = None
                effective_type = typ
                if section.kind == task_node.SectionType.RUN:
                    effective_type, parsed_default = _normalize_run_param(name, typ, parsed.default)
                param = Param(name, effective_type, parsed.default, parsed.desc, parsed_default)
                if section.kind == task_node.SectionType.IN:
                    inputs[param.name] = param
                elif section.kind == task_node.SectionType.OUT:
                    if param.default is None:
                        raise ValueError(f'Output parameter must have a default: {param.name}')
                    outputs[param.name] = param
                elif section.kind == task_node.SectionType.RUN:
                    run[param.name] = param
                else:
                    raise ValueError(f'Unrecognized section kind: {section.kind}')

    return TaskSignature(inputs, outputs, run)


def _normalize_run_param(name: str, typ: TypeKind, default):
    entry = {
        'memory': (TypeKind.MEMORY, 'memory', _parse_memory_literal),
        'time': (TypeKind.TIME, 'time', _parse_time_literal),
        'cpu': (TypeKind.INT, 'int', _parse_cpu_literal),
        'image': (TypeKind.STR, 'str', None),
    }.get(name)
    if entry is None:
        return typ, None
    expected_kind, expected_label, parser = entry
    if typ is None:
        typ = expected_kind
    elif typ != expected_kind:
        raise ValueError(f'Run parameter {name} must have type {expected_label}: {typ.value}')
    if default is None:
        return typ, None
    text = _literal_default_text(default)
    if text is not None:
        value = text if parser is None else parser(text)
        return typ, value
    if _has_interpolation(default):
        return typ, None
    raise ValueError(f'Run parameter {name} must have a literal default')


def _literal_default_text(default):
    if not isinstance(default, interp.Word):
        return None
    if len(default.parts) != 1:
        return None
    part = default.parts[0]
    if not isinstance(part, interp.Literal):
        return None
    return part.text.strip()


def _has_interpolation(default):
    if not isinstance(default, interp.Word):
        return False
    return any(isinstance(part, interp.Var) for part in default.parts)


def _parse_cpu_literal(text: str) -> int:
    if not re.fullmatch(r'[0-9]+', text):
        raise ValueError(f'Invalid cpu literal: {text}')
    return int(text)


def _parse_memory_literal(text: str) -> int:
    m = re.fullmatch(r'([0-9]+)([KMGTP])?', text)
    if not m:
        raise ValueError(f'Invalid memory literal: {text}')
    value = int(m.group(1))
    unit = m.group(2) or 'M'
    scale = {
        'K': 1 / 1024,
        'M': 1,
        'G': 1024,
        'T': 1024 * 1024,
        'P': 1024 * 1024 * 1024,
    }[unit]
    return int(math.ceil(value * scale))


def _parse_time_literal(text: str) -> int:
    if re.fullmatch(r'[0-9]+', text):
        return int(text)
    m = re.fullmatch(r'(?:(\d+)-)?(\d+):(\d+):(\d+)', text)
    if not m:
        raise ValueError(f'Invalid time literal: {text}')
    days = int(m.group(1) or 0)
    hours = int(m.group(2))
    minutes = int(m.group(3))
    seconds = int(m.group(4))
    total_seconds = days * 24 * 60 * 60 + hours * 60 * 60 + minutes * 60 + seconds
    return int(math.ceil(total_seconds / 60))
