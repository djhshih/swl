import traceback

from swl.ir.lower import lower_file
from swl.ir import node as ir


def _print(node, indent=''):
    if isinstance(node, ir.Block):
        print(f'{indent}Block')
        for bind in node.bindings:
            print(f'{indent}  Bind {bind.name}:')
            _print(bind.value, indent + '    ')
        print(f'{indent}  Result:')
        _print(node.result, indent + '    ')
        return

    if isinstance(node, ir.Bind):
        print(f'{indent}Bind {node.name}')
        _print(node.value, indent + '  ')
        return

    if isinstance(node, ir.Function):
        print(f'{indent}Function name={node.name} kind={node.kind} path={node.path}')
        print(f'{indent}  inputs={sorted(node.signature.inputs.keys())}')
        print(f'{indent}  outputs={sorted(node.signature.outputs.keys())}')
        if node.body is not None:
            print(f'{indent}  body:')
            _print(node.body, indent + '    ')
        return

    if isinstance(node, ir.Lambda):
        print(f'{indent}Lambda param={node.param}')
        _print(node.body, indent + '  ')
        return

    if isinstance(node, ir.Record):
        print(f'{indent}Record')
        for name, value in node.fields.items():
            print(f'{indent}  {name}:')
            _print(value, indent + '    ')
        return

    if isinstance(node, ir.Field):
        print(f'{indent}Field {node.name}')
        _print(node.record, indent + '  ')
        return

    if isinstance(node, ir.Update):
        print(f'{indent}Update')
        _print(node.left, indent + '  ')
        _print(node.right, indent + '  ')
        return

    if isinstance(node, ir.Apply):
        print(f'{indent}Apply')
        print(f'{indent}  function:')
        _print(node.function, indent + '    ')
        print(f'{indent}  arg:')
        _print(node.arg, indent + '    ')
        return

    if isinstance(node, ir.Chain):
        print(f'{indent}Chain')
        for item in node.items:
            _print(item, indent + '  ')
        return

    if isinstance(node, ir.Name):
        print(f'{indent}Name {node.name}')
        return

    if isinstance(node, ir.Literal):
        print(f'{indent}Literal {node.value!r}')
        return

    if isinstance(node, ir.Unknown):
        print(f'{indent}Unknown')
        return

    print(f'{indent}{node!r}')


def eval(fname):
    tree = lower_file(fname)
    _print(tree)


if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Lower workflow to semantic IR')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)
