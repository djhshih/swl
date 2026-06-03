import argparse
import json


def run(description, ext, emit_fn, is_json):
    ap = argparse.ArgumentParser(description)
    ap.add_argument('input', help='compiled DAG json path')
    ap.add_argument('-o', '--output', help=f'output .{ext} path')
    args = ap.parse_args()
    result = emit_fn(args.input)
    payload = json.dumps(result, indent=2) if is_json else result
    if args.output:
        with open(args.output, 'w') as f:
            f.write(payload)
            f.write('\n')
    else:
        print(payload)
