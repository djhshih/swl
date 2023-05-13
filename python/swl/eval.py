import traceback

import swl.syntax.parser as pr
import swl.syntax.lexer as lr
from swl.syntax.node import NodeType

def print_ast(node, depth=0):
    # fun and block nodes need to be indented
    indent = ' ' * depth
    if node.type == NodeType.fun:
        print(f'{indent}(\\ {node.param}')
        print_ast(node.body, depth + 1)
        print(f'{indent})')
    elif node.type == NodeType.block:
        print(f'{indent}[')
        for child in node.body:
            print_ast(child, depth + 1)
        print(f'{indent}]')
    else:
        print(f'{indent}{node}')

def eval(fname):
    p = pr.Parser()
    with open(fname, 'r') as f:
        src = f.read()
        tokens = [t for t in lr.Lexer(src)]
        print('tokens:')
        print(tokens)
        tree = p.parse(src)
        print('tree:')
        print_ast(tree)
    

if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Parse simple workflow language')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)

