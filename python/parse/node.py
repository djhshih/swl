from enum import Enum

NodeType = Enum('NodeType',
    [
        'id',
        'str',
        'num',
        'block',
        'bind',
        'rec',
        'fun',
        'apply',
        'get',
        'pipe',
    ]
)

class Expr:
    type = None

class Binding(Expr):
    type = NodeType.bind
    def __init__(self, id, expr):
        self.id = id
        self.value = expr

class Block(Expr):
    type = NodeType.block
    def __init__(self, lines):
        self.lines = lines

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

class Get(Expr):
    type = NodeType.get
    def __init__(self, rec, member):
        self.rec = rec
        self.member = member

class Pipe(Expr):
    type = NodeType.pipe
    def __init__(self, funs):
        self.funs = funs

class Record(Expr):
    type = NodeType.rec
    def __init__(self, dict):
        self.dict = dict

class String(Expr):
    type = NodeType.str
    def __init__(self, value):
        self.value = value

class Number(Expr):
    type = NodeType.num
    def __init__(self, value):
        self.value = value

class Identifier(Expr):
    type = NodeType.id
    def __init__(self, name):
        self.name = name

