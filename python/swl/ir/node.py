from dataclasses import dataclass, field
from typing import Dict, List, Optional

from swl.semantic.task.type import TaskSignature


@dataclass(frozen=True)
class Node:
    pass


@dataclass(frozen=True)
class Literal(Node):
    value: object


@dataclass(frozen=True)
class Unknown(Node):
    pass


@dataclass(frozen=True)
class Name(Node):
    name: str


@dataclass(frozen=True)
class Record(Node):
    fields: Dict[str, Node]


@dataclass(frozen=True)
class Field(Node):
    record: Node
    name: str


@dataclass(frozen=True)
class Update(Node):
    left: Node
    right: Node


@dataclass(frozen=True)
class Function(Node):
    name: str
    kind: str
    signature: TaskSignature
    path: Optional[str] = None
    body: Optional[Node] = None


@dataclass(frozen=True)
class Lambda(Node):
    param: str
    body: Node
    signature: Optional[TaskSignature] = None


@dataclass(frozen=True)
class Closure(Node):
    function: Node
    bound_arg: Node
    signature: Optional[TaskSignature] = None


@dataclass(frozen=True)
class Apply(Node):
    function: Node
    arg: Node
    signature: Optional[TaskSignature] = None


@dataclass(frozen=True)
class Chain(Node):
    items: List[Node]
    signature: Optional[TaskSignature] = None


@dataclass(frozen=True)
class Stage:
    name: str
    function: Node
    arg: Node


@dataclass(frozen=True)
class Compose(Node):
    param: Optional[str]
    stages: List[Stage]
    result: Node
    signature: Optional[TaskSignature] = None


@dataclass(frozen=True)
class Bind(Node):
    name: str
    value: Node


@dataclass(frozen=True)
class Block(Node):
    bindings: List[Bind] = field(default_factory=list)
    result: Node = field(default_factory=Unknown)
