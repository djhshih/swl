import argparse
import json

from swl.transpile.wdl.emit import transpile_dag_file


def main():
    ap = argparse.ArgumentParser('Transpile compiled DAG JSON to WDL 1.1')
    ap.add_argument('input', help='compiled DAG json path')
    ap.add_argument('-o', '--output', help='output .wdl path')
    args = ap.parse_args()

    wdl = transpile_dag_file(args.input)
    if args.output:
        with open(args.output, 'w') as f:
            f.write(wdl)
            f.write('\n')
    else:
        print(wdl)


if __name__ == '__main__':
    main()
