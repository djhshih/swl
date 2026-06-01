import argparse
import json

from swl.transpile.nf.emit import transpile_dag_file


def main():
    ap = argparse.ArgumentParser('Transpile compiled DAG JSON to Nextflow DSL2')
    ap.add_argument('input', help='compiled DAG json path')
    ap.add_argument('-o', '--output', help='output .nf path')
    args = ap.parse_args()

    nf = transpile_dag_file(args.input)
    if args.output:
        with open(args.output, 'w') as f:
            f.write(nf)
            f.write('\n')
    else:
        print(nf)


if __name__ == '__main__':
    main()
