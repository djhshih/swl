import argparse
import os
import sys
import traceback

from swl.api import compile_workflow


if __name__ == '__main__':
    ap = argparse.ArgumentParser('Compile workflow to canonical logical DAG JSON')
    ap.add_argument('input', help='input workflow file')
    ap.add_argument('-o', '--output', help='output DAG path')
    args = ap.parse_args()

    try:
        output = compile_workflow(args.input, args.output)
        print(f'wrote canonical logical DAG: {output}')
    except Exception:
        traceback.print_exc(file=sys.stderr)
        raise SystemExit(os.EX_DATAERR)
