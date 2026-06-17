# SWL Testing Upgrade Plan

## 1. Trivial & Restrictive Tests

### Trivial Tests (insufficient depth)

These tests verify the bare minimum and miss edge cases:

| File | Issue |
|------|-------|
| `tests/unit/swl/syntax/task/test_interpolation.py` | Only 5 tests, each testing one trivial happy-path. Missing: empty input, nested `${}`, special chars, multiple vars, escaped `$`, invalid syntax variations |
| `tests/unit/swl/syntax/task/test_bash.py` | Only 2 tests parsing 2 script bodies. Missing: empty scripts, comments in body, heredocs, pipes, redirects, error recovery |
| `tests/unit/swl/syntax/wf/test_lexer.py` | 7 tests, all happy-path token sequences. Missing: error cases (unterminated strings, unexpected chars, empty input, deeply nested `{}`) |
| `tests/unit/swl/syntax/wf/test_parser.py` | 10 of 14 tests use `assertIsNotNone` — they check parsing doesn't crash but never inspect the AST structure emitted. Missing: AST shape verification, error recovery, complex nested structures |
| `tests/unit/swl/semantic/wf/test_validate.py` | Only 4 tests covering 2 workflows. Missing: missing required inputs, wrong types, extra unknown inputs, optional inputs, mixed batch/simple |
| `tests/unit/swl/transpile/nf/test_emit.py` | `test_empty_inputs_outputs`, `test_outputspec_passthrough` — single assertion per test |
| `tests/unit/swl/transpile/wdl/test_emit.py` | `test_empty_inputs_outputs` — trivial smoke test |

### Restrictive Tests (brittle, impede refactoring)

These assert exact structure/strings that make the code hard to refactor:

| File | Test | Issue |
|------|------|-------|
| `test_force.py` | `test_serialized_dag_is_self_contained` | Asserts exact path `data['steps'][0]['path']` |
| `test_force.py` | `test_force_saturated_workflow_produces_task_dag` | Hardcodes step IDs `['align', 'sort']` and output keys |
| `test_force.py` | `test_mapped_table_source_uses_explicit_logical_table_metadata` | Asserts entire 5-column dict verbatim |
| `test_force.py` | `test_chain_root_is_instantiated_during_force` | 16 assertions on exact structure |
| `test_lower.py` | `test_lower_imports_to_functions` | 13 assertions on exact IR node structure chain |
| `test_lower.py` | `test_chain_and_explicit_have_equivalent_lowered_shape` | Good in spirit but fragile — any IR shape change breaks both |
| CWL `test_emit.py` | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` | Asserts 10 scatter port names verbatim |
| CWL `test_emit.py` | `test_root_partial_map_transpiles_as_scattered_subworkflow` | Asserts 10 scatter port names verbatim |
| NF `test_emit.py` | `test_root_partial_map_inlines_pipeline` | Asserts exact string `'SORT(ALIGN.out.bam, outbase)'` |
| NF `test_emit.py` | `test_process_name_sanitization` | Asserts `_process_name('_align') == 'ALIGN'` — losing `_` prefix silently |

---

## 2. Coverage Gaps — Untested Modules

### Critical — Entire modules with zero tests:

| Module | Functions/Methods | Risk |
|--------|-------------------|------|
| `swl/api.py` | `compile_workflow`, `force_workflow`, `load_workflow`, `transpile_dag` | Public API is completely untested |
| `swl/compile.py` | CLI arg parsing, error handling | Entry-point has zero coverage |
| `swl/loader.py` | `Loader.__init__`, `read_file`, `get_parsed_task`, `cache_task`, `get_checked_workflow`, `cache_workflow` | Core file I/O caching untested |
| `swl/repl.py` | REPL loop | REPL is untested |
| `swl/transpile/smk/emit.py` | All 25+ functions (`transpile_dag_dict`, `_task_to_rule`, `_interpolate_shell`, `_dag_to_smk`, etc.) | **Complete Snakemake backend** has zero unit tests |
| `swl/transpile/common.py` | `normalize_identifier`, `workflow_name`, `step_name`, `emit_name`, `field_chain_parts`, `field_path_after_first`, `classify_var`, `interp_script`, `word_interp`, `run_value`, `source_kind`, etc. | Shared utilities used by ALL backends — untested |
| `swl/transpile/_cli.py` | `run()` | CLI runner untested |
| `swl/dag/merge.py` | `_canonicalize_merges`, `_flatten_value_terms`, `_normalize_output_value`, `_value_key`, `_merge_key`, `_record_fields_key` | Merge normalization logic untested |
| `swl/dag/tooldefs.py` | All 15+ functions (`_tool_definition`, `_workflow_definition`, `_materialize_workflow_dag`, `_force_root`, etc.) | Tool/workflow definition logic untested |
| `swl/eval/ir.py` | `eval()` | IR evaluator untested |
| `swl/eval/dag.py` | `eval()` | DAG evaluator untested |
| `swl/eval/semantic_wf.py` | `eval()`, `_format_param` | Workflow semantic evaluator untested |
| `swl/eval/semantic_task.py` | `eval()`, `_format_param` | Task semantic evaluator untested |
| `swl/eval/syntax_task.py` | `eval()` | Task syntax evaluator untested |
| `swl/eval/syntax_wf.py` | `eval()`, `print_ast` | Workflow syntax evaluator untested |

### All transpile CLI/`__main__` modules (zero tests):
- `swl/transpile/cwl/cli.py`, `cwl/__main__.py`
- `swl/transpile/nf/cli.py`, `nf/__main__.py`
- `swl/transpile/wdl/cli.py`, `wdl/__main__.py`
- `swl/transpile/smk/cli.py`, `smk/__main__.py`

---

## 3. Coverage Gaps — Untested Functions Within Tested Modules

### `swl/ir/lower.py`
- `_lower_inline_import`
- `_generated_callable_from_lambda`
- Other private `_lower_*` helpers

### `swl/dag/forcer.py`
- `_is_opaque_record_carrier`
- `_ensure_forced_record`
- `make_force_state` constructor logic (tested indirectly but never directly)

### `swl/dag/evaluator.py`
- Most of `_force_map`, `_apply`, `_available_inputs`, `_value_key`, `force_value` internals (only tested end-to-end)

### `swl/dag/node.py`
- `DAG.validate()` — only error paths tested (circular, self-loop, unknown dep); success path not tested
- `DAG.from_dict()` — only a few `from_dict` round-trip tests; many node types not validated
- `Field`, `Record`, `Literal`, `Merge` serialization/deserialization not tested

### `swl/dag/context.py`
- `ForceEnv` — all methods untested

### `swl/dag/binding.py`
- All binding utilities untested

### `swl/dag/finalize.py`
- `_prune_unused_inputs`, `_build_output_specs`, `_finalize_dag` — only tested indirectly

### `swl/semantic/wf/infer.py`
- `application_result`, `apply_function`, `ClosureValue`, `ComputationValue`, `ClosedRecord` — most inference logic only indirectly tested

### `swl/semantic/wf/bashvars.py`
- `_validate_bash_variables` — only indirectly tested when it throws; no direct tests

### `swl/semantic/wf/scope.py`
- Scope checking entirely untested

### `swl/semantic/wf/imports.py`
- Import resolution logic untested

### `swl/semantic/wf/signature.py`
- Signature operations untested

### `swl/syntax/wf/builtins.py`
- `map`, `map_by` builtin logic untested

### `swl/semantic/wf/type.py`
- `FunctionType`, `RecordType`, `TableType` — partially tested through check tests but type constructors/operations not directly tested

---

## 4. Recommended Test Plan (Priority Ordered)

### Phase 1 — Shore up existing tests (low effort, high impact)

1. **Deepen interpolation tests** (`test_interpolation.py`): Add edge cases — empty `${}`, nested `${}`, escaped `\$`, multiple vars, invalid syntax
2. **Deepen lexer tests** (`test_lexer.py`): Add error recovery — unterminated strings, unexpected chars, empty string, max nesting
3. **Add AST structure checks** to `test_parser.py`: Don't just `assertIsNotNone` — verify node types, children, positions
4. **Deepen bash parser tests** (`test_bash.py`): Add heredocs, redirects, pipelines, comments, empty body, error handling
5. **Deepen validation tests** (`test_validate.py`): Missing inputs, wrong types, extra inputs, optional handling

### Phase 2 — Add unit tests for shared utilities (medium effort, high impact)

1. **`swl/transpile/common.py`**: Unit test all 12+ functions. Especially `interp_script`, `word_interp`, `normalize_identifier`, `field_chain_parts`
2. **`swl/dag/merge.py`**: Unit test `_canonicalize_merges`, `_normalize_output_value`, `_value_key` with various merge tree shapes
3. **`swl/dag/node.py`**: Unit test `DAG.validate()` success path, all `from_dict` node types, `StepCall` construction

### Phase 3 — Add Snakemake backend tests (medium effort, high impact)

1. **Create `tests/unit/swl/transpile/smk/test_emit.py`**: Test `transpile_dag_dict`, `_task_to_rule`, `_dag_to_smk`, `_interpolate_shell`, `_binding_to_path`, `_collect_params`, `_emit_resources`, `_validate_supported`
2. Match coverage pattern of existing CWL/NF/WDL test suites

### Phase 4 — API and Loader tests (medium effort, high impact)

1. **`swl/api.py`**: Unit test `compile_workflow`, `force_workflow`, `load_workflow`, `transpile_dag` with in-memory virtual filesystem
2. **`swl/loader.py`**: Unit test caching behavior, read_file from both virtual files and real filesystem

### Phase 5 — Core evaluation and forcing internals (high effort, high impact)

1. **`swl/dag/evaluator.py`**: Direct unit tests for `_force_map`, `_apply`, `_available_inputs` with controlled IR inputs
2. **`swl/dag/tooldefs.py`**: Unit test `_tool_definition`, `_workflow_definition`, `_force_root`, `_materialize_workflow_dag`
3. **`swl/dag/forcer.py`**: Direct test for `_is_opaque_record_carrier`, `make_force_state`

### Phase 6 — Semantic inference and analysis (medium effort, medium impact)

1. **`swl/semantic/wf/infer.py`**: Direct unit tests for `apply_function`, `application_result` with mock signatures
2. **`swl/semantic/wf/scope.py`**: Test scope resolution, shadow rules, forward references
3. **`swl/semantic/wf/imports.py`**: Test circular import detection, caching, path resolution
4. **`swl/semantic/wf/bashvars.py`**: Test `_validate_bash_variables` with known/unknown variables

### Phase 7 — Lowerer internals (low effort, medium impact)

1. **`swl/ir/lower.py`**: Direct tests for `_lower_inline_import`, `_generated_callable_from_lambda`
2. **`swl/syntax/wf/builtins.py`**: Direct tests for `map`/`map_by` builtin resolution

### Phase 8 — REPL and CLIs (low effort, low impact)

1. **`swl/repl.py`**: Basic smoke test with simulated stdin
2. **All `cli.py`/`__main__.py`**: Test argument parsing and error handling

---

## 5. Summary Metrics

| Category | Count |
|----------|-------|
| Python source modules | 33 (excluding empty `__init__.py`) |
| Modules with any unit tests | ~14 |
| Modules with zero tests | ~19 |
| Functions with zero tests | ~100+ |
| Trivial/brittle tests needing hardening | ~20 |
