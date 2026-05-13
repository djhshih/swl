import traceback

import swl.syntax.task.parser as pr
import swl.semantic.task.type as ty


def _format_param(param):
    extra = ''
    if param.parsed_default is not None:
        extra = f' parsed_default={param.parsed_default!r}'
    return f'type={param.type} default={param.default!r}{extra} desc={param.desc!r}'


def eval(fname):
    p = pr.Parser()
    with open(fname, 'r') as f:
        src = f.read()
        task = p.parse(src)
        sig = ty.signature_from_task(task)
        print('annotation:')
        print(f'doc: {task.annotation.doc}')
        print('signature:')
        print('inputs:')
        for name, param in sig.inputs.items():
            print(f'  {name}: {_format_param(param)}')
        print('outputs:')
        for name, param in sig.outputs.items():
            print(f'  {name}: {_format_param(param)}')
        print('run:')
        for name, param in sig.run.items():
            print(f'  {name}: {_format_param(param)}')


if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Parse task annotation and build semantic signature')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)
