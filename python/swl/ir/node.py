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
class Ref(Node):
    id: int
    name: str

    def __repr__(self):
        return f'Ref(id={self.id!r}, name={self.name!r})'


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

    def __repr__(self):
        return (
            f'Field(\n'
            f'  record={self.record!r},\n'
            f'  name={self.name!r},\n'
            f')'
        )


@dataclass(frozen=True)
class ArrayField(Node):
    record_array: Node
    name: str

    def __repr__(self):
        return (
            f'ArrayField(\n'
            f'  record_array={self.record_array!r},\n'
            f'  name={self.name!r},\n'
            f')'
        )


@dataclass(frozen=True)
class Update(Node):
    left: Node
    right: Node

    def __repr__(self):
        return (
            f'Update(\n'
            f'  left={self.left!r},\n'
            f'  right={self.right!r},\n'
            f')'
        )


@dataclass(frozen=True)
class Function(Node):
    name: str
    kind: str
    signature: TaskSignature
    path: Optional[str] = None
    body: Optional[Node] = None

    def __repr__(self):
        return (
            f'Function(\n'
            f'  name={self.name!r},\n'
            f'  kind={self.kind!r},\n'
            f'  signature={self.signature!r},\n'
            f'  path={self.path!r},\n'
            f'  body={self.body!r},\n'
            f')'
        )


@dataclass(frozen=True)
class Lambda(Node):
    param: str
    body: 'Block'
    signature: Optional[TaskSignature] = None

    def __repr__(self):
        return (
            f'Lambda(\n'
            f'  param={self.param!r},\n'
            f'  body={self.body!r},\n'
            f'  signature={self.signature!r},\n'
            f')'
        )


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

    def __repr__(self):
        return (
            f'Apply(\n'
            f'  function={self.function!r},\n'
            f'  arg={self.arg!r},\n'
            f'  signature={self.signature!r},\n'
            f')'
        )


@dataclass(frozen=True)
class Map(Node):
    function: Node
    arg: Node

    def __repr__(self):
        return (
            f'Map(\n'
            f'  function={self.function!r},\n'
            f'  arg={self.arg!r},\n'
            f')'
        )


@dataclass(frozen=True)
class Variable(Node):
    id: int
    name: str
    value: Node

    def __repr__(self):
        return (
            f'Variable(\n'
            f'  id={self.id!r},\n'
            f'  name={self.name!r},\n'
            f'  value={self.value!r},\n'
            f')'
        )


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
