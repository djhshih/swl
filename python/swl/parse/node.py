from enum import Enum

NodeType = Enum('NodeType',
    [
        'id',
        'str',
        'num',
        'block',
        'bind',
        'rec',
        'get',
        'update',
        'fun',
        'apply',
        'pipe',
    ]
)

class Expr:
    type = None
    def __repr__(self):
        if self.type == NodeType.id:
            return f'{self.name}'
        elif self.type == NodeType.str:
            return f'"{self.value}"'
        elif self.type == NodeType.num:
            return f'{self.value}'
        elif self.type == NodeType.block:
            return '[' + ';\n'.join([str(x) for x in self.body]) + ']'
        elif self.type == NodeType.bind:
            return f'(= {self.id} {self.value})'
        elif self.type == NodeType.rec:
            return ','.join([f'{k}: {v}' for k, v in self.value])
        elif self.type == NodeType.get:
            return f'(. {self.rec} {self.member})'
        elif self.type == NodeType.update:
            return f'(& {self.left} {self.right})'
        elif self.type == NodeType.fun:
            return f'(\ {self.param} {self.body})'
        elif self.type == NodeType.apply:
            return f'($ {self.fun} {self.arg})'
        elif self.type == NodeType.pipe:
            return f'(|> {self.left} {self.right})'
        else:
            return f'{self.type}'

class Identifier(Expr):
    type = NodeType.id
    def __init__(self, name):
        self.name = name

class String(Expr):
    type = NodeType.str
    def __init__(self, value):
        self.value = value

class Number(Expr):
    type = NodeType.num
    def __init__(self, value):
        self.value = value

class Block(Expr):
    type = NodeType.block
    def __init__(self, lines):
        self.body = lines

class Binding(Expr):
    type = NodeType.bind
    def __init__(self, iden, expr):
        self.id = iden
        self.value = expr

class Record(Expr):
    type = NodeType.rec
    def __init__(self, dict):
        self.value = dict

class Get(Expr):
    type = NodeType.get
    def __init__(self, rec, member):
        self.rec = rec
        self.member = member

class Update(Expr):
    type = NodeType.update
    def __init__(self, left, right):
        self.left = left
        self.right = right

class Function(Expr):
    type = NodeType.fun
    def __init__(self, param, block):
        self.param = param
        self.body = block

class Apply(Expr):
    type = NodeType.apply
    def __init__(self, fun, arg):
        self.fun = fun
        self.arg = arg

class Pipe(Expr):
    type = NodeType.pipe
    def __init__(self, left, right):
        self.left = left
        self.right = right

