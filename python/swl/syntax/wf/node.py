from enum import Enum

NodeType = Enum('NodeType', [
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
    'chain',
])

REPR_FORMATS = {
    NodeType.id: lambda self: str(getattr(self, 'name', '<missing>')),
    NodeType.str: lambda self: str(getattr(self, 'value', '<missing>')),
    NodeType.num: lambda self: str(getattr(self, 'value', '<missing>')),
    NodeType.block: lambda self: '[' + '; '.join([f'{x}' for x in getattr(self, 'body', [])]) + ']',
    NodeType.bind: lambda self: f'(= {getattr(self, "id", "<missing>")} {getattr(self, "value", "<missing>")})',
    NodeType.rec: lambda self: '{' + ', '.join([f'{k}: {v}' for k, v in getattr(self, 'value', {}).items()]) + '}',
    NodeType.get: lambda self: f'(. {getattr(self, "rec", "<missing>")} {getattr(self, "member", "<missing>")})',
    NodeType.update: lambda self: f'(// {getattr(self, "left", "<missing>")} {getattr(self, "right", "<missing>")})',
    NodeType.fun: lambda self: f'(\\ {getattr(self, "param", "<missing>")} {getattr(self, "body", "<missing>")})',
    NodeType.apply: lambda self: f'($ {getattr(self, "fun", "<missing>")} {getattr(self, "arg", "<missing>")})',
    NodeType.chain: lambda self: f'(| {getattr(self, "left", "<missing>")} {getattr(self, "right", "<missing>")})',
}

class Expr:
    type = None
    def __repr__(self):
        fmt = REPR_FORMATS.get(self.type)
        if fmt is not None:
            return fmt(self)
        return f'{self.type}'

class Identifier(Expr):
    type = NodeType.id
    def __init__(self, name):
        self.name = name
        self.binding = None

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
    def __init__(self, d):
        self.value = d

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
    def __init__(self, param, body):
        self.param = param
        self.body = body

class Apply(Expr):
    type = NodeType.apply
    def __init__(self, fun, arg):
        self.fun = fun
        self.arg = arg

class Chain(Expr):
    type = NodeType.chain
    def __init__(self, left, right):
        self.left = left
        self.right = right

