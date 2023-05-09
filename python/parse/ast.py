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

class ExprN:
    type = None

class BindingN:
    type = NodeType.bind
    def __init__(self, id, expr):
        self.id = id
        self.value = expr

class BlockN:
    type = NodeType.block
    def __init__(self, lines):
        self.lines = lines

class FunctionN:
    type = NodeType.fun
    def __init__(self, param, block):
        self.param = param
        self.body = block

class ApplyN:
    type = NodeType.apply
    def __init__(self, fun, arg):
        self.fun = fun
        self.arg = arg

class GetN:
    type = NodeType.get
    def __init__(self, rec, member):
        self.rec = rec
        self.member = member

class PipeN:
    type = NodeType.pipe
    def __init__(self, funs):
        self.funs = funs

class RecordN:
    type = NodeType.rec
    def __init__(self, dict):
        self.dict = dict

class StringN:
    type = NodeType.str
    def __init__(self, value):
        self.value = value

class NumberN:
    type = NodeType.num
    def __init__(self, value):
        self.value = value

class IdentifierN:
    type = NodeType.id
    def __init__(self, name):
        self.name = name

