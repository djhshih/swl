from swl.syntax.task import interpolation as interp
from swl.syntax.task import bash as task_bash


def _bash_interp_vars(word):
    vars = []
    parts = word.parts if isinstance(word, interp.Word) else [word]
    for part in parts:
        if isinstance(part, interp.Var):
            vars.append(part.name)
    return vars


def _validate_bash_variables(parsed_body, input_names, context_name):
    errors = []
    defined = set(input_names)
    for stmt in parsed_body.statements:
        if isinstance(stmt, task_bash.Assignment):
            refs = _bash_interp_vars(stmt.value)
            for ref in refs:
                if ref not in defined:
                    errors.append(
                        f'Unresolved variable "${ref}" in assignment to "{stmt.name}" '
                        f'in {context_name}'
                    )
            defined.add(stmt.name)
        elif isinstance(stmt, task_bash.Command):
            for word in stmt.words:
                for ref in _bash_interp_vars(word):
                    if ref not in defined:
                        errors.append(
                            f'Unresolved variable "${ref}" in command "{stmt.text[:60]}" '
                            f'in {context_name}'
                        )
    return errors
