import argparse
import os
import sys
import traceback

from swl.dag.forcer import force_file


def compile_workflow(input_path, output_path=None):
    dag = force_file(input_path)
    if output_path is None:
        base = os.path.splitext(os.path.basename(input_path))[0]
        output_path = os.path.join(os.path.dirname(input_path), 'dag', f'{base}.json')
    dag.write(output_path)
    return output_path


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
