import math
import re
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


def parse_type(type_str: str) -> TypeKind:
    '''Parse a type string into a TypeKind.'''
    type_str = type_str.strip()

    type_map = {
        'file': TypeKind.FILE,
        'str': TypeKind.STR,
        'int': TypeKind.INT,
        'float': TypeKind.FLOAT,
        'memory': TypeKind.MEMORY,
        'time': TypeKind.TIME,
        'file?': TypeKind.FILE_OPT,
        'str?': TypeKind.STR_OPT,
        'int?': TypeKind.INT_OPT,
        'float?': TypeKind.FLOAT_OPT,
        '[file]': TypeKind.ARRAY_FILE,
        '[str]': TypeKind.ARRAY_STR,
        '[int]': TypeKind.ARRAY_INT,
        '[float]': TypeKind.ARRAY_FLOAT,
    }

    if type_str not in type_map:
        raise ValueError(f'Unknown type: {type_str}')

    return type_map[type_str]


def is_optional(typ: TypeKind) -> bool:
    return typ in (
        TypeKind.FILE_OPT,
        TypeKind.STR_OPT,
        TypeKind.INT_OPT,
        TypeKind.FLOAT_OPT,
    )


def is_array(typ: TypeKind) -> bool:
    return typ in (
        TypeKind.ARRAY_FILE,
        TypeKind.ARRAY_STR,
        TypeKind.ARRAY_INT,
        TypeKind.ARRAY_FLOAT,
    )


def base_type(typ: TypeKind) -> TypeKind:
    base_map = {
        TypeKind.FILE_OPT: TypeKind.FILE,
        TypeKind.STR_OPT: TypeKind.STR,
        TypeKind.INT_OPT: TypeKind.INT,
        TypeKind.FLOAT_OPT: TypeKind.FLOAT,
        TypeKind.ARRAY_FILE: TypeKind.FILE,
        TypeKind.ARRAY_STR: TypeKind.STR,
        TypeKind.ARRAY_INT: TypeKind.INT,
        TypeKind.ARRAY_FLOAT: TypeKind.FLOAT,
    }
    return base_map.get(typ, typ)


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


class Param:
    '''Semantic task parameter.'''

    def __init__(self, name: str, typ: TypeKind | None = None, default: str | None = None, desc: str | None = None, parsed_default=None):
        self.name = name
        self.type = typ
        self.default = default
        self.desc = desc
        self.parsed_default = parsed_default


class TaskSignature:
    '''Type signature for a task.'''

    def __init__(self, inputs: Dict[str, Param], outputs: Dict[str, Param], run: Dict[str, Param] = None):
        self.inputs = inputs
        self.outputs = outputs
        self.run = run or {}

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

    def infer_workflow_input(self, tasks: List[str]) -> Set[str]:
        all_outputs = set()
        for task in tasks:
            if task in self.signatures:
                all_outputs.update(self.signatures[task].output_names())

        unsatisfied = set()
        for task in tasks:
            if task in self.signatures:
                inputs = self.signatures[task].input_names()
                for inp in inputs:
                    if inp not in all_outputs:
                        unsatisfied.add(inp)

        return unsatisfied


def signature_from_task(task: task_node.Task) -> TaskSignature:
    '''Build semantic task signature from parsed task syntax.'''
    inputs = {}
    outputs = {}
    run = {}

    for section in task.annotation.sections:
        for parsed in section.params:
            typ = parse_type(parsed.type) if parsed.type else None
            for name in parsed.names:
                parsed_default = None
                effective_type = typ
                if section.kind == task_node.SectionType.RUN:
                    effective_type, parsed_default = _normalize_run_param(name, typ, parsed.default)
                param = Param(name, effective_type, parsed.default, parsed.desc, parsed_default)
                if section.kind == task_node.SectionType.IN:
                    _add_param(inputs, param, 'input')
                elif section.kind == task_node.SectionType.OUT:
                    if param.default is None:
                        raise ValueError(f'Output parameter must have a default: {param.name}')
                    _add_param(outputs, param, 'output')
                elif section.kind == task_node.SectionType.RUN:
                    _add_param(run, param, 'run')
                else:
                    raise ValueError(f'Unrecognized section kind: {section.kind}')

    return TaskSignature(inputs, outputs, run)


def _add_param(params: Dict[str, Param], param: Param, kind: str):
    if param.name in params:
        raise ValueError(f'Duplicate {kind} parameter: {param.name}')
    params[param.name] = param


def _normalize_run_param(name: str, typ: TypeKind, default):
    if name == 'memory':
        return _normalize_memory_param(typ, default)
    if name == 'time':
        return _normalize_time_param(typ, default)
    if name == 'cpu':
        return _normalize_cpu_param(typ, default)
    if name == 'image':
        return _normalize_image_param(typ, default)
    return typ, None


def _normalize_memory_param(typ: TypeKind, default):
    if typ is None:
        typ = TypeKind.MEMORY
    elif typ != TypeKind.MEMORY:
        raise ValueError(f'Run parameter memory must have type memory: {typ.value}')
    if default is None:
        return typ, None
    text = _literal_default_text(default)
    if text is None:
        raise ValueError('Run parameter memory must have a literal default')
    return typ, _parse_memory_literal(text)


def _normalize_time_param(typ: TypeKind, default):
    if typ is None:
        typ = TypeKind.TIME
    elif typ != TypeKind.TIME:
        raise ValueError(f'Run parameter time must have type time: {typ.value}')
    if default is None:
        return typ, None
    text = _literal_default_text(default)
    if text is None:
        raise ValueError('Run parameter time must have a literal default')
    return typ, _parse_time_literal(text)


def _normalize_cpu_param(typ: TypeKind, default):
    if typ is None:
        typ = TypeKind.INT
    elif typ != TypeKind.INT:
        raise ValueError(f'Run parameter cpu must have type int: {typ.value}')
    if default is None:
        return typ, None
    text = _literal_default_text(default)
    if text is None:
        raise ValueError('Run parameter cpu must have a literal default')
    return typ, _parse_cpu_literal(text)


def _normalize_image_param(typ: TypeKind, default):
    if typ is None:
        typ = TypeKind.STR
    elif typ != TypeKind.STR:
        raise ValueError(f'Run parameter image must have type str: {typ.value}')
    if default is None:
        return typ, None
    text = _literal_default_text(default)
    if text is None:
        raise ValueError('Run parameter image must have a literal default')
    return typ, text


def _literal_default_text(default):
    if not isinstance(default, interp.Word):
        return None
    if len(default.parts) != 1:
        return None
    part = default.parts[0]
    if not isinstance(part, interp.Literal):
        return None
    return part.text.strip()


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
