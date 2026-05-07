import traceback

import swl.semantic.wf.check as ck


def eval(fname):
    result = ck.Checker().load(fname)
    print('imports:')
    for name, imported in result.imports.items():
        print(f'  {name}: kind={imported.kind} path={imported.path}')
    print('chain errors:')
    print(result.chain_errors)
    print('inferred inputs:')
    print(sorted(result.inferred_inputs))
    print('signature:')
    if result.signature is None:
        print('  None')
    else:
        print('  inputs:')
        for name, param in result.signature.inputs.items():
            print(f'    {name}: type={param.type} default={param.default!r} desc={param.desc!r}')
        print('  outputs:')
        for name, param in result.signature.outputs.items():
            print(f'    {name}: type={param.type} default={param.default!r} desc={param.desc!r}')
        print('  run:')
        for name, param in result.signature.run.items():
            print(f'    {name}: type={param.type} default={param.default!r} desc={param.desc!r}')


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
