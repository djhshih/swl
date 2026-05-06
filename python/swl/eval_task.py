import traceback

import swl.syntax.task.parser as pr


def eval(fname):
    p = pr.Parser()
    with open(fname, 'r') as f:
        src = f.read()
        task = p.parse(src)
        print('annotation:')
        print(f'doc: {task.annotation.doc}')
        for section in task.annotation.sections:
            print(f'section: {section.kind.value}')
            for param in section.params:
                print(
                    '  param:',
                    {
                        'names': param.names,
                        'type': param.type,
                        'default': repr(param.default),
                        'desc': param.desc,
                    }
                )
        print('body:')
        print(task.body)


if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Parse task annotation and body')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)
