import json
import traceback

from swl.dag.forcer import force_file


def eval(fname):
    dag = force_file(fname)
    print('canonical logical DAG:')
    print(f'  inputs: {sorted(dag.inputs.keys())}')
    print(f'  steps: {[step.id for step in dag.steps]}')
    mapped = [step for step in dag.steps if getattr(step, 'map', None) is not None]
    if mapped:
        print('  mapped steps:')
        for step in mapped:
            source = getattr(step, 'map', {}).get('source', {})
            source_kind = source.get('source')
            if source_kind == 'table':
                source_desc = f"table(columns={sorted(source.get('columns', {}).keys())})"
            elif source_kind == 'input':
                source_desc = f"input({source.get('name')})"
            else:
                source_desc = source_kind or 'unknown'
            print(f'    {step.id}: type={step.type} source={source_desc}')
            print(f'      input_schema={getattr(step, "input_schema", None)}')
            print(f'      output_schema={getattr(step, "output_schema", None)}')
    print(json.dumps(dag.to_dict(), indent=2, sort_keys=True))


if __name__ == '__main__':
    import sys, os, argparse

    ap = argparse.ArgumentParser('Force workflow to canonical logical DAG')
    ap.add_argument('input', help='input file')

    args = ap.parse_args()

    try:
        eval(args.input)
    except:
        traceback.print_exc()
        sys.exit(os.EX_DATAERR)
