import traceback

from swl.ir.lower import Lowerer
import swl.ir.node as ir


def eval(fname):
    lowerer = Lowerer()
    tree = lowerer.lower_file(fname)
    print('semantic IR:')
    print('functions:')
    for name in sorted(lowerer.function_cache.keys()):
        for line in repr(lowerer.function_cache[name]).splitlines():
            print(f'  {line}')
    print('tree:')
    if isinstance(tree, ir.Function):
        print(f'  Function {tree.name}')
    else:
        for line in repr(tree).splitlines():
            print(f'  {line}')


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
