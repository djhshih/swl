import os

from swl.semantic.task.type import Param, TaskSignature, signature_from_task
from swl.semantic.wf.bashvars import _validate_bash_variables
from swl.syntax.task import bash as task_bash
from swl.syntax.task.parser import Parser as TaskParser
from swl.syntax.wf import builtins, node as wf_node


class Import:
    def __init__(self, name: str, path: str, signature: TaskSignature, kind: str, check=None, task=None, parsed_body=None):
        self.name = name
        self.path = path
        self.signature = signature
        self.kind = kind
        self.check = check
        self.task = task
        self.parsed_body = parsed_body


def read_file(checker, path: str) -> str:
    if path in checker.files:
        return checker.files[path]
    with open(path, 'r') as f:
        return f.read()


def load_imports(checker, tree, base_dir: str):
    imports = {}
    if tree.type != wf_node.NodeType.block:
        return imports
    for expr in tree.body:
        if expr.type != wf_node.NodeType.bind:
            continue
        path = builtins.match_import(expr.value)
        if path is None:
            continue
        full_path = os.path.join(base_dir, path)
        imports[expr.id.name] = load_import(checker, expr.id.name, full_path)
    return imports


def load_import(checker, name: str, path: str) -> Import:
    if path.endswith('.sh'):
        src = read_file(checker, path)
        task = TaskParser().parse(src)
        parsed_body = task_bash.Parser().parse(task.body)
        signature = signature_from_task(task)
        known_vars = set(signature.inputs.keys()) | set(signature.run.keys())
        var_errors = _validate_bash_variables(parsed_body, known_vars, f'task "{name}" ({path})')
        if var_errors:
            raise ValueError('\n'.join(var_errors))
        return Import(name, path, signature, 'task', task=task, parsed_body=parsed_body)
    if path.endswith('.swl'):
        check = checker.load(path)
        if check.signature is None:
            raise ValueError(f'Imported workflow does not produce a signature: {path}')
        return Import(name, path, check.signature, 'workflow', check=check)
    raise ValueError(f'Unrecognized import path: {path}')
