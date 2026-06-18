from swl.semantic.scope import DuplicateBindingError, Scope
from swl.syntax.wf import builtins, node as wf_node


def check_scope(checker, tree):
    errors = []
    walk_scope(checker, tree, Scope(), errors)
    return errors


def walk_scope(checker, expr, outer_scope, errors):
    if expr.type == wf_node.NodeType.block:
        block_scope = Scope(parent=outer_scope)
        bind_order = []
        for item in expr.body:
            if item.type == wf_node.NodeType.bind:
                bind_order.append(item)
                if builtins.match_import(item.value) is None:
                    walk_scope(checker, item.value, block_scope, errors)
                try:
                    block_scope.declare(item.id.name, node=item)
                except DuplicateBindingError:
                    errors.append(f'Duplicate binding in scope: {item.id.name}')
            else:
                walk_scope(checker, item, block_scope, errors)
        for i, item in enumerate(bind_order):
            defined_so_far = _collect_names(outer_scope)
            for j in range(i):
                defined_so_far.add(bind_order[j].id.name)
            if builtins.match_import(item.value) is None:
                refs = collect_name_refs(checker, item.value)
            else:
                refs = []
            own_name = item.id.name
            for ref in refs:
                if ref == own_name:
                    errors.append(
                        f'Forward reference in binding scope: {own_name} '
                        f'references itself'
                    )
                elif ref not in defined_so_far:
                    is_future = any(
                        ref == later.id.name
                        for later in bind_order[i+1:]
                    )
                    if is_future:
                        errors.append(
                            f'Forward reference in binding scope: {own_name} '
                            f'references {ref} before its definition'
                        )
        return

    if expr.type == wf_node.NodeType.fun:
        fn_scope = Scope(parent=outer_scope)
        fn_scope.declare(expr.param.name)
        walk_scope(checker, expr.body, fn_scope, errors)
        return

    for child in children(checker, expr):
        walk_scope(checker, child, outer_scope, errors)


def _collect_names(scope):
    names = set(scope.locals.keys())
    while scope.parent is not None:
        scope = scope.parent
        names.update(scope.locals.keys())
    return names


def collect_name_refs(checker, expr):
    refs = []
    walk_refs(checker, expr, refs)
    return refs


def walk_refs(checker, expr, refs):
    if expr.type == wf_node.NodeType.id:
        refs.append(expr.name)
        return
    if expr.type == wf_node.NodeType.fun:
        return
    for child in children(checker, expr):
        walk_refs(checker, child, refs)


def check_chains(checker, tree, type_checker):
    errors = []
    walk_chains(checker, tree, type_checker, errors)
    return errors


def walk_chains(checker, expr, type_checker, errors):
    if expr.type == wf_node.NodeType.chain:
        left_name = task_name(checker, expr.left)
        right_name = task_name(checker, expr.right)
        if left_name and right_name:
            errors.extend(type_checker.check_chain(left_name, right_name))
    for child in children(checker, expr):
        walk_chains(checker, child, type_checker, errors)


def chain_names(checker, expr):
    if expr.type == wf_node.NodeType.id:
        return [expr.name]
    if expr.type == wf_node.NodeType.chain:
        return chain_names(checker, expr.left) + chain_names(checker, expr.right)
    return []


def task_name(checker, expr):
    if expr.type == wf_node.NodeType.id:
        return expr.name
    if expr.type == wf_node.NodeType.chain:
        names = chain_names(checker, expr)
        if names:
            return names[-1]
    return None


def annotate_ast(checker, tree):
    _walk_annotate(checker, tree, Scope())


def _walk_annotate(checker, expr, scope):
    if expr.type == wf_node.NodeType.id:
        binding = scope.resolve(expr.name)
        if binding is not None:
            expr.binding = binding
        return

    if expr.type == wf_node.NodeType.block:
        block_scope = Scope(parent=scope)
        for item in expr.body:
            if item.type == wf_node.NodeType.bind:
                if builtins.match_import(item.value) is None:
                    _walk_annotate(checker, item.value, block_scope)
                block_scope.declare(item.id.name, node=item)
            else:
                _walk_annotate(checker, item, block_scope)
        return

    if expr.type == wf_node.NodeType.fun:
        fn_scope = Scope(parent=scope)
        fn_scope.declare(expr.param.name)
        _walk_annotate(checker, expr.body, fn_scope)
        return

    for child in children(checker, expr):
        _walk_annotate(checker, child, scope)


def children(checker, expr):
    if expr.type == wf_node.NodeType.block:
        return expr.body
    if expr.type == wf_node.NodeType.bind:
        return [expr.id, expr.value]
    if expr.type == wf_node.NodeType.rec:
        return list(expr.value.values())
    if expr.type == wf_node.NodeType.get:
        return [expr.rec, expr.member]
    if expr.type == wf_node.NodeType.update:
        return [expr.left, expr.right]
    if expr.type == wf_node.NodeType.fun:
        return [expr.param, expr.body]
    if expr.type == wf_node.NodeType.apply:
        return [expr.fun, expr.arg]
    if expr.type == wf_node.NodeType.chain:
        return [expr.left, expr.right]
    return []
