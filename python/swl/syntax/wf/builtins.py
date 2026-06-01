from swl.syntax.wf import node as wf_node


def match_import(expr):
    if expr.type != wf_node.NodeType.apply:
        return None
    if expr.fun.type != wf_node.NodeType.id or expr.fun.name != 'import':
        return None
    if expr.arg.type != wf_node.NodeType.str:
        return None
    return expr.arg.value


def match_map(expr):
    if expr.type != wf_node.NodeType.apply:
        return None
    left = expr.fun
    if left.type != wf_node.NodeType.apply:
        return None
    if left.fun.type != wf_node.NodeType.id or left.fun.name != 'map':
        return None
    return left.arg, expr.arg


def match_map_by(expr):
    if expr.type != wf_node.NodeType.apply:
        return None
    left = expr.fun
    if left.type != wf_node.NodeType.apply:
        return None
    inner = left.fun
    if inner.type != wf_node.NodeType.apply:
        return None
    if inner.fun.type != wf_node.NodeType.id or inner.fun.name != 'map_by':
        return None
    return left.arg, inner.arg, expr.arg
