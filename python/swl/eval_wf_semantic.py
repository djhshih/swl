import traceback

import swl.semantic.wf.check as ck


def _format_param(param):
    extra = ''
    if getattr(param, 'parsed_default', None) is not None:
        extra = f' parsed_default={param.parsed_default!r}'
    return f'type={param.type} default={param.default!r}{extra} desc={param.desc!r}'


def eval(fname):
    result = ck.Checker().load(fname)
    print('imports:')
    for name, imported in result.imports.items():
        print(f'  {name}: kind={imported.kind} path={imported.path}')
    print('errors:')
    print(result.errors)
    print('inferred inputs:')
    print(sorted(result.inferred_inputs))
    print('signature:')
    if result.signature is None:
        print('  None')
    else:
        print('  inputs:')
        for name, param in result.signature.inputs.items():
            print(f'    {name}: {_format_param(param)}')
        print('  outputs:')
        for name, param in result.signature.outputs.items():
            print(f'    {name}: {_format_param(param)}')
        print('  run:')
        for name, param in result.signature.run.items():
            print(f'    {name}: {_format_param(param)}')
    print(f'workflow_type: {result.workflow_type!r}')
    print(f'root_input_type: {result.root_input_type!r}')
    print(f'root_output_type: {result.root_output_type!r}')


if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Parse workflow and build semantic import view')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)
