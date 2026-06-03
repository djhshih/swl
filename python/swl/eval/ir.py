import traceback

from swl.ir.lower import Lowerer
from swl.ir import node as ir


def _print(node, indent='', tree_mode=False):
    if isinstance(node, ir.Function) and tree_mode:
        print(f'{indent}Function {node.name}')
        return

    text = repr(node)
    for line in text.splitlines():
        print(f'{indent}{line}')


def eval(fname):
    lowerer = Lowerer()
    tree = lowerer.lower_file(fname)
    print('semantic IR:')
    print('functions:')
    for name in sorted(lowerer.function_cache.keys()):
        _print(lowerer.function_cache[name], '  ', False)
    print('tree:')
    _print(tree, '  ', True)


if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Lower workflow to semantic IR (pre-force, pre-CWL)')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)
