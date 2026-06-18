from dataclasses import dataclass, field, fields
from typing import Dict, List, Optional, Set

from swl.semantic.task.type import TaskSignature


@dataclass(frozen=True)
class Node:
    _repr_fields = None

    def __repr__(self):
        names = self._repr_fields
        if names is None:
            names = tuple(f.name for f in fields(self))
        if not names:
            return f'{type(self).__name__}()'
        parts = ''.join(f'  {name}={getattr(self, name)!r},\n' for name in names)
        return f'{type(self).__name__}(\n{parts})'


@dataclass(frozen=True)
class Literal(Node):
    value: object


@dataclass(frozen=True)
class Input(Node):
    name: str
    type: Optional[str] = None
    desc: Optional[str] = None


@dataclass(frozen=True)
class Unknown(Node):
    pass


@dataclass(frozen=True)
class Name(Node):
    name: str


@dataclass(frozen=True)
class Ref(Node):
    id: int
    name: str


@dataclass(frozen=True)
class Record(Node):
    fields: Dict[str, Node]

    def __repr__(self):
        if self.fields:
            fields = ',\n'.join(
                f'    {name!r}: {value!r}'
                for name, value in self.fields.items()
            )
            fields = f'{{\n{fields}\n  }}'
        else:
            fields = '{}'
        return f'Record(fields={fields})'


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
    is_batch: bool = False


@dataclass(frozen=True)
class Lambda(Node):
    param: str
    body: 'Block'
    signature: Optional[TaskSignature] = None
    is_batch: bool = False
    captures: Dict[str, Node] = field(default_factory=dict)


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
    satisfied: Set[str] = field(default_factory=set)

    _repr_fields = ('function', 'arg', 'signature')


@dataclass(frozen=True)
class Map(Node):
    function: Node
    arg: Node
    key: Optional[str] = None


@dataclass(frozen=True)
class Variable(Node):
    id: int
    name: str
    value: Node


@dataclass(frozen=True)
class Block(Node):
    bindings: List[Variable] = field(default_factory=list)
    result: Node = field(default_factory=Unknown)

    def __repr__(self):
        if self.bindings:
            bindings = ',\n'.join(f'    {bind!r}' for bind in self.bindings)
            bindings = f'[\n{bindings}\n  ]'
        else:
            bindings = '[]'
        return (
            f'Block(\n'
            f'  bindings={bindings},\n'
            f'  result={self.result!r},\n'
            f')'
        )
