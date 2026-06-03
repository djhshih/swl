from swl.ir import node as ir
from swl.ir.lower import Lowerer
from swl.dag.binding import binding_to_dict
from swl.dag.context import ForceEnv
from swl.dag.evaluator import (
    _assert_canonical,
    _available_inputs,
    _is_opaque_record_carrier,
    _looks_record_like,
    _can_project_record_fields,
    _project_field,
    _project_input,
    _merge_values,
    _merge_bound,
    _register_block,
    _saturated_signature,
    _saturated_task,
    _forced_apply_key,
    _apply_key,
    force_value,
    _force_ref,
    _force_apply,
    _apply,
    _force_lambda,
    _force_map,
    _apply_mapped,
    _is_batch_function,
)
from swl.dag.emit import (
    _emit_task_call,
    _emit_workflow_call,
    _emit_mapped_step,
    _mapped_step_bindings,
    _normalize_task_inputs,
    _normalize_task_run,
    _step_id,
    _step_call_key,
    _validate_bindings,
    _step_dependencies,
    _walk_values,
    _value_dependencies,
    _logical_table_source,
    _binding_to_public_dict,
)
from swl.dag.finalize import (
    _finalize_dag,
    _refine_input_metadata,
    _flatten_step_bindings,
    _flatten_outputs,
    _build_output_specs,
    _prune_unused_inputs,
    _assert_wireable_output,
    _infer_output_type,
    _final_outputs,
    _collect_output_fields,
    _flatten_merge_value,
    _unwrap_non_record,
    _merge_input_type,
    _merge_input_desc,
)
from swl.dag.tooldefs import (
    _tool_definition,
    _workflow_definition,
    _materialize_workflow_dag,
    _build_task_section,
    _build_task_param,
    _validate_output_default_glob,
    _input,
    _forced_signature,
    _function_signature,
    _force_root,
    _is_unsaturated_builtin_map,
    _materialize_partial_map_workflow,
)
from swl.types import is_optional_type, normalize_swl_type, to_array_type


class Forcer:
    def __init__(self, files=None):
        self.inputs = {}
        self.steps = []
        self.step_cache = {}
        self.tool_defs = {}
        self.call_counter = 0
        self.step_name_counts = {}
        self.lowerer = Lowerer(files=files)
        self.variables = {}
        self.forced_variables = {}
        self.apply_cache = {}

    def force(self, node):
        _assert_canonical(self, node)
        value = force_value(self, node, ForceEnv())
        value = _force_root(self, value)
        return _finalize_dag(self, value)

    def _finalize_dag(self, value):
        return _finalize_dag(self, value)

    def _assert_canonical(self, node):
        _assert_canonical(self, node)

    def force_value(self, node, env):
        return force_value(self, node, env)

    def _force_ref(self, ref, env):
        return _force_ref(self, ref, env)

    def _force_apply(self, node, env):
        return _force_apply(self, node, env)

    def _apply(self, fn, arg, satisfied=None):
        return _apply(self, fn, arg, satisfied)

    def _force_lambda(self, function, bound):
        return _force_lambda(self, function, bound)

    def _force_map(self, fn, arg, key=None):
        return _force_map(self, fn, arg, key=key)

    def _apply_mapped(self, fn, source, key=None):
        return _apply_mapped(self, fn, source, key=key)

    def _is_batch_function(self, fn):
        return _is_batch_function(self, fn)

    def _register_block(self, bindings):
        _register_block(self, bindings)

    def _merge_bound(self, old, new):
        return _merge_bound(self, old, new)

    def _available_inputs(self, value):
        return _available_inputs(self, value)

    def _saturated_signature(self, function, bound):
        return _saturated_signature(self, function, bound)

    def _saturated_task(self, function, bound):
        return _saturated_task(self, function, bound)

    def _forced_apply_key(self, value):
        return _forced_apply_key(self, value)

    def _apply_key(self, function_node, fn, arg):
        return _apply_key(self, function_node, fn, arg)

    def _looks_record_like(self, value):
        return _looks_record_like(self, value)

    def _can_project_record_fields(self, value):
        return _can_project_record_fields(self, value)

    def _is_opaque_record_carrier(self, value):
        return _is_opaque_record_carrier(self, value)

    def _project_field(self, source, name):
        return _project_field(self, source, name)

    def _merge_values(self, left, right):
        return _merge_values(self, left, right)

    def _project_input(self, value, name):
        return _project_input(self, value, name)

    def _normalize_task_inputs(self, function, bound):
        return _normalize_task_inputs(self, function, bound)

    def _emit_task_call(self, function, bound, satisfied=None):
        return _emit_task_call(self, function, bound, satisfied)

    def _emit_workflow_call(self, function, bound, satisfied=None):
        return _emit_workflow_call(self, function, bound, satisfied)

    def _emit_mapped_step(self, fn, source, key=None):
        return _emit_mapped_step(self, fn, source, key=key)

    def _mapped_step_bindings(self, fn, source, key=None):
        return _mapped_step_bindings(self, fn, source, key=key)

    def _step_id(self, name):
        return _step_id(self, name)

    def _step_call_key(self, function, inputs):
        return _step_call_key(self, function, inputs)

    def _validate_bindings(self, inputs, satisfied, function_name):
        return _validate_bindings(self, inputs, satisfied, function_name)

    def _step_dependencies(self, inputs):
        return _step_dependencies(self, inputs)

    def _walk_values(self, value):
        return _walk_values(self, value)

    def _value_dependencies(self, value):
        return _value_dependencies(self, value)

    def _logical_table_source(self, source):
        return _logical_table_source(self, source)

    def _binding_to_public_dict(self, value):
        return _binding_to_public_dict(value)

    def _refine_input_metadata(self):
        _refine_input_metadata(self)

    def _flatten_step_bindings(self):
        _flatten_step_bindings(self)

    def _flatten_outputs(self, outputs):
        return _flatten_outputs(self, outputs)

    def _build_output_specs(self, outputs):
        return _build_output_specs(self, outputs)

    def _prune_unused_inputs(self, root_value, outputs):
        _prune_unused_inputs(self, root_value, outputs)

    def _assert_wireable_output(self, name, value):
        return _assert_wireable_output(self, name, value)

    def _infer_output_type(self, value):
        return _infer_output_type(self, value)

    def _final_outputs(self, value):
        return _final_outputs(self, value)

    def _collect_output_fields(self, value):
        return _collect_output_fields(self, value)

    def _flatten_merge_value(self, value):
        return _flatten_merge_value(self, value)

    def _unwrap_non_record(self, value):
        return _unwrap_non_record(self, value)

    def _merge_input_type(self, name, current, candidates):
        return _merge_input_type(self, name, current, candidates)

    def _merge_input_desc(self, current, candidates):
        return _merge_input_desc(self, current, candidates)

    def _tool_definition(self, path):
        return _tool_definition(self, path)

    def _workflow_definition(self, function):
        return _workflow_definition(self, function)

    def _materialize_workflow_dag(self, function):
        return _materialize_workflow_dag(self, function)

    def _build_task_section(self, task, signature, kind):
        return _build_task_section(self, task, signature, kind)

    def _build_task_param(self, signature, kind, name, param):
        return _build_task_param(self, signature, kind, name, param)

    def _validate_output_default_glob(self, name, default):
        return _validate_output_default_glob(self, name, default)

    def _input(self, name, spec=None):
        return _input(self, name, spec)

    def _forced_signature(self, value):
        return _forced_signature(self, value)

    def _function_signature(self, function):
        return _function_signature(self, function)

    def _force_root(self, value):
        return _force_root(self, value)

    def _is_unsaturated_builtin_map(self, value, signature):
        return _is_unsaturated_builtin_map(self, value, signature)

    def _materialize_partial_map_workflow(self, value, signature):
        return _materialize_partial_map_workflow(self, value, signature)

    def _normalize_task_run(self, function, bound):
        return _normalize_task_run(self, function, bound)


def force_file(path: str, files=None):
    forcer = Forcer(files=files)
    tree = forcer.lowerer.lower_file(path)
    return forcer.force(tree)
