from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass(frozen=True)
class Input:
    name: str
    type: Optional[str] = None
    desc: Optional[str] = None


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
    signature: Optional[object] = None


@dataclass
class DAG:
    inputs: Dict[str, Input]
    tasks: List[TaskCall]
    outputs: Dict[str, object]
