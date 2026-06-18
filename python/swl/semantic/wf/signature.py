from swl.semantic.scope import Scope
from swl.semantic.task.type import Param, TaskSignature
from swl.semantic.wf import type as wf_type
from swl.semantic.wf.infer import (
    ClosedRecord,
    ClosureValue,
    ComputationValue,
    FunctionValue,
    OpenRecord,
    TableValue,
    UnknownValue,
    eval_function_body,
    eval_lambda_in_mode,
    eval_prefix_bindings,
    infer_root_record_type,
    infer_root_table_type,
    initial_param_value,
    inner_partial_remaining_inputs,
    is_function_value,
    looks_like_stale_record_errors,
    signature_inputs_from_root_type,
    signature_inputs_from_table_type,
    signature_outputs,
    value_is_table,
    value_type,
)
from swl.syntax.wf import builtins, node as wf_node


def build_workflow_signature(checker, tree, imports, inferred_inputs, issues):
    if tree.type != wf_node.NodeType.block or not tree.body:
        return None, False, None, None, None

    final = tree.body[-1]
    if final.type == wf_node.NodeType.fun:
        env = eval_prefix_bindings(checker, tree.body[:-1], imports, set(), issues)

        from swl.semantic.wf.infer import uses_map
        if uses_map(checker, final.body):
            demanded, mode_issues, chosen_result = eval_lambda_in_mode(checker, final, imports, env, 'table', inferred_inputs)
            if value_is_table(chosen_result):
                issues[:] = [
                    issue for issue in issues
                    if not issue.startswith('Cannot resolve field access: (. ')
                    and issue not in ('map requires a tab argument', 'map_by requires a tab argument')
                ]
                root_input_type = infer_root_table_type(checker, final, imports, env, inferred_inputs)
                root_output_type = wf_output_type(checker, chosen_result, table_context=True)
                inputs = signature_inputs_from_table_type(checker, root_input_type)
                if not inputs:
                    inputs = {final.param.name: Param(final.param.name, None)}
                outputs = signature_outputs(checker, chosen_result)
                signature = TaskSignature(inputs, outputs, {})
                return signature, True, wf_type.FunctionType(root_input_type, root_output_type), root_input_type, root_output_type
            filtered_issues = [] if looks_like_stale_record_errors(mode_issues) else mode_issues
            issues.extend(filtered_issues)
            root_input_type = infer_root_table_type(checker, final, imports, env, inferred_inputs)
            root_output_type = wf_output_type(checker, chosen_result, table_context=True)
            inputs = signature_inputs_from_table_type(checker, root_input_type)
            if not inputs:
                inputs = {final.param.name: Param(final.param.name, None)}
            outputs = signature_outputs(checker, chosen_result)
            signature = TaskSignature(inputs, outputs, {})
            return signature, True, wf_type.FunctionType(root_input_type, root_output_type), root_input_type, root_output_type

        record_scope = Scope(parent=env)
        record_scope.declare(final.param.name).value = initial_param_value('record', inferred_inputs)
        chosen_result = eval_function_body(checker, final, imports, record_scope, set(), list(issues))
        root_output_type = wf_output_type(checker, chosen_result, table_context=False)
        partial_inputs = inner_partial_remaining_inputs(checker, final, imports, env)
        if isinstance(chosen_result, ClosureValue):
            inputs = dict(chosen_result.signature.inputs)
            root_input_type = wf_type.RecordType({
                name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
                for name, param in sorted(inputs.items())
            }, open=True)
        elif partial_inputs is not None:
            inputs = dict(partial_inputs)
            root_input_type = wf_type.RecordType({
                name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
                for name, param in sorted(inputs.items())
            }, open=True)
        else:
            root_input_type = infer_root_record_type(checker, final, imports, env, inferred_inputs)
            inputs = signature_inputs_from_root_type(checker, root_input_type)
        outputs = signature_outputs(checker, chosen_result)
        signature = TaskSignature(inputs, outputs, {})
        return signature, False, wf_type.FunctionType(root_input_type, root_output_type), root_input_type, root_output_type

    scope = Scope()
    for expr in tree.body[:-1]:
        if expr.type == wf_node.NodeType.bind:
            if builtins.match_import(expr.value) is not None:
                continue
            from swl.semantic.wf.infer import eval_expr
            scope.set_local(expr.id.name, value=eval_expr(checker, expr.value, imports, scope, set(), issues))
    from swl.semantic.wf.infer import eval_expr
    result = eval_expr(checker, final, imports, scope, set(), issues)
    if isinstance(result, FunctionValue):
        wft = wf_function_type_from_signature(checker, result.signature, result.batch)
        return result.signature, result.batch, wft, getattr(wft, 'input', None), getattr(wft, 'output', None)
    if isinstance(result, ClosureValue):
        signature = result.signature
        wft = wf_function_type_from_signature(checker, signature, result.function.batch)
        return signature, result.function.batch, wft, getattr(wft, 'input', None), getattr(wft, 'output', None)
    return None, False, None, None, None


def wf_output_type(checker, value, table_context=False):
    if isinstance(value, TableValue):
        return wf_type.TableType({name: wf_type.UNKNOWN for name in sorted(value.columns.keys())})
    if isinstance(value, (ClosedRecord, OpenRecord)):
        fields = getattr(value, 'fields', {})
        return wf_type.RecordType({name: wf_type.UNKNOWN for name in sorted(fields.keys())}, open=isinstance(value, OpenRecord))
    if isinstance(value, (ComputationValue, ClosureValue, FunctionValue)):
        outputs = signature_outputs(checker, value)
        return wf_type.RecordType({name: wf_type.UNKNOWN for name in sorted(outputs.keys())}, open=False)
    return wf_type.RecordType({}, open=True)


def wf_function_type_from_signature(checker, signature, is_batch):
    if signature is None:
        return None
    output = wf_type.RecordType({
        name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
        for name, param in sorted(signature.outputs.items())
    }, open=False)
    if is_batch:
        input_type = wf_type.TableType({
            name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
            for name, param in sorted(signature.inputs.items())
        })
    else:
        input_type = wf_type.RecordType({
            name: wf_type.scalar_from_name(getattr(getattr(param, 'type', None), 'value', None))
            for name, param in sorted(signature.inputs.items())
        }, open=True)
    return wf_type.FunctionType(input_type, output)
