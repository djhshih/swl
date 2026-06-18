import os

from swl.semantic.task.type import Param, TaskSignature, signature_from_task
from swl.semantic.wf.bashvars import _validate_bash_variables
from swl.syntax.task import bash as task_bash, interpolation as interp
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
    return checker.loader.read_file(path)


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
        default_errors = _validate_param_defaults(signature.run, known_vars, f'task "{name}" ({path})')
        default_errors.extend(_validate_param_defaults(signature.outputs, known_vars, f'task "{name}" ({path})'))
        var_errors.extend(default_errors)
        if var_errors:
            raise ValueError('\n'.join(var_errors))
        checker.loader.cache_task(path, task, signature, parsed_body)
        return Import(name, path, signature, 'task', task=task, parsed_body=parsed_body)
    if path.endswith('.swl'):
        check = checker.load(path)
        if check.signature is None:
            raise ValueError(f'Imported workflow does not produce a signature: {path}')
        return Import(name, path, check.signature, 'workflow', check=check)
    raise ValueError(f'Unrecognized import path: {path}')


def _validate_param_defaults(params, known_vars, context_name):
    errors = []
    for name, param in params.items():
        if param.default is None:
            continue
        default_word = param.default
        if isinstance(default_word, interp.Word):
            for part in default_word.parts:
                if isinstance(part, interp.Var):
                    if part.name not in known_vars:
                        errors.append(
                            f'Unresolved variable "${{{part.name}}}" in {name} default '
                            f'in {context_name}'
                        )
                elif isinstance(part, interp.Expr):
                    for var_name in task_bash._extract_expr_vars(part.text):
                        if var_name not in known_vars:
                            errors.append(
                                f'Unresolved variable "${var_name}" in {name} default '
                                f'in {context_name}'
                            )
    return errors
