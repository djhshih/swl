from swl.syntax.wf import node as wf_node


def check_scope(checker, tree):
    errors = []
    walk_scope(checker, tree, set(), errors)
    return errors


def walk_scope(checker, expr, scope, errors):
    if expr.type == wf_node.NodeType.block:
        local = set(scope)
        bind_order = []
        for item in expr.body:
            if item.type == wf_node.NodeType.bind:
                bind_order.append(item)
                if item.id.name in local:
                    errors.append(f'Duplicate binding in scope: {item.id.name}')
                else:
                    local.add(item.id.name)
            else:
                walk_scope(checker, item, local, errors)
        for i, item in enumerate(bind_order):
            defined_so_far = set(scope)
            for j in range(i):
                defined_so_far.add(bind_order[j].id.name)
            refs = collect_name_refs(checker, item.value)
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
            walk_scope(checker, item.value, local, errors)
        return

    if expr.type == wf_node.NodeType.fun:
        local = {expr.param.name}
        walk_scope(checker, expr.body, local, errors)
        return

    for child in children(checker, expr):
        walk_scope(checker, child, scope, errors)


def collect_name_refs(checker, expr):
    refs = []
    walk_refs(checker, expr, refs)
    return refs


def walk_refs(checker, expr, refs):
    if expr.type == wf_node.NodeType.id:
        refs.append(expr.name)
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
