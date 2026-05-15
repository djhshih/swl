import argparse
import json

from swl.transpile.cwl.emit import transpile_dag_file


def main():
    ap = argparse.ArgumentParser('Transpile compiled DAG JSON to packed CWL')
    ap.add_argument('input', help='compiled DAG json path')
    ap.add_argument('-o', '--output', help='output CWL path')
    args = ap.parse_args()

    data = transpile_dag_file(args.input)
    payload = json.dumps(data, indent=2)
    if args.output:
        with open(args.output, 'w') as f:
            f.write(payload)
            f.write('\n')
    else:
        print(payload)


if __name__ == '__main__':
    main()
