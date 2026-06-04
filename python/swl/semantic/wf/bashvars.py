from swl.syntax.task import bash as task_bash


def _validate_bash_variables(parsed_body, known_vars, context_name):
    errors = []
    defined = set(known_vars)
    for stmt in parsed_body.statements:
        if isinstance(stmt, task_bash.Assignment):
            for name, _ in task_bash.iter_var_refs(task_bash.Script([stmt])):
                if name not in defined and name not in task_bash._BUILTIN_VARS:
                    errors.append(
                        f'Unresolved variable "${name}" in assignment to "{stmt.name}" '
                        f'in {context_name}'
                    )
            defined.add(stmt.name)
        elif isinstance(stmt, task_bash.Command):
            for name, _ in task_bash.iter_var_refs(task_bash.Script([stmt])):
                if name not in defined and name not in task_bash._BUILTIN_VARS:
                    errors.append(
                        f'Unresolved variable "${name}" in command "{stmt.text[:60]}" '
                        f'in {context_name}'
                    )
    return errors
