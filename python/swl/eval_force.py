import json
import traceback

from swl.ir.force import force_file


def eval(fname):
    dag = force_file(fname)
    print(json.dumps(dag.to_dict(), indent=2, sort_keys=True))


if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Force workflow to execution DAG')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)
