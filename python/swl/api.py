import json
import os

from swl.dag.forcer import force_file as _force_file, make_force_state
from swl.ir.lower import Lowerer as _Lowerer


def compile_workflow(input_path, output_path=None):
    dag = _force_file(input_path)
    if output_path is None:
        base = os.path.splitext(os.path.basename(input_path))[0]
        output_path = os.path.join(os.path.dirname(input_path), 'dag', f'{base}.json')
    dag.write(output_path)
    return output_path


def force_workflow(path, files=None):
    state = make_force_state(files=files)
    tree = state.lowerer.lower_file(path)
    from swl.dag.forcer import force
    return force(tree, state=state)


def load_workflow(path, files=None):
    lowerer = _Lowerer(files=files)
    return lowerer.checker.load(path)


def transpile_dag(dag_path, target):
    target = target.lower()
    if target == 'cwl':
        from swl.transpile.cwl.emit import transpile_dag_file
        return json.dumps(transpile_dag_file(dag_path), indent=2)
    if target == 'wdl':
        from swl.transpile.wdl.emit import transpile_dag_file
        return transpile_dag_file(dag_path)
    if target == 'nf':
        from swl.transpile.nf.emit import transpile_dag_file
        return transpile_dag_file(dag_path)
    if target == 'smk':
        from swl.transpile.smk.emit import transpile_dag_file
        return transpile_dag_file(dag_path)
    raise ValueError(f'Unsupported transpilation target: {target!r}')
