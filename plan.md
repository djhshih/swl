# Remaining implementation plan

## Goal

Drive `python/swl/` toward production readiness. Tracks concrete implementation work derived from:
- `dag.md` — the DAG JSON spec and cross-transpiler requirements
- `cwl.md` — CWL transpiler status and gaps
- `wdl.md` — WDL 1.1 transpiler implementation plan
- `nf.md` — Nextflow transpiler implementation plan
- `spec.md` — SWL language specification
- `progress.md` — current gap list

---



## Priority 1: CWL transpiler — remaining gaps

Driven by `cwl.md`. Gaps that remain after `map_by`, `expr` passthrough, and time hints are implemented.

### 1a. Record bindings via ExpressionTool

**Status:** ✅ Complete. Record bindings in step inputs and workflow outputs are handled by `_emit_record_tool` (emit.py:519) which generates an `ExpressionTool`. Tests:
- `test_record_binding_in_step_input_emits_expression_tool` (test_emit.py:378)
- `test_record_workflow_output_emits_expression_tool` (test_emit.py:436)
- `test_record_binding_no_longer_rejected` (test_emit.py:474)

### 1b. Nested `Field(Field(...))` projections

**Status:** ✅ Complete. The CWL transpiler handles nested `Field(Field(...))` chains via `_canonical_binding` in emit.py:307-329, which produces `input_field_nested`, `step_output_nested`, and `tab_column_step_output_nested` kinds. These get a `valueFrom` expression in `_step_input_to_cwl`. Test:
- `test_nested_field_binding_uses_valueFrom` (test_emit.py:496)

### 1c. Merge workflow outputs

**Status:** ✅ Complete. Merge flattening at the output level is fully handled:
- `_canonicalize_merges` in force.py:969 collapses record-record merges
- `_flatten_outputs` in force.py:661 calls `_flatten_merge_value` which recursively flattens merge trees
- `_unwrap_non_record` now preserves `Field(source=StepCall)` values (was silently dropping them)
- End-to-end tests exist via `test_function_and_chain_compile_to_same_shape` (test_force.py:325) which uses `//` at output level
- New P1c-specific tests added:
  - `test_merge_record_field_stepcall_preserves_field` — Field(StepCall) no longer dropped
  - `test_merge_record_field_input_preserves_field` — Field(Input) preserved
  - `test_merge_record_with_bare_input_drops_input` — bare Input still dropped (intentional)
  - `test_end_to_end_merge_at_output_level_leaves_no_merge_in_dag` — no Merge in final DAG

---

## Priority 2: Table update semantics

**Driven by:** `spec.md §Record Update`, `progress.md`

### 2a. Table-record updates (`tab // rec`, `rec // tab`)

**Status:** The semantic checker (`_merge_table_record` in `check.py`) handles table-record scalar broadcast. The force phase raises `ValueError` when a `StepCall` appears on either side of `//`. The DAG construction path for table field-level updates during forcing is not implemented.

**What to do:**
- Review the spec requirement: "If a table `t` is updated with a record `r` via `t // r`, scalar properties in `r` are implicitly duplicated across all rows of `t`."
- Identify whether the gap is in the spec, the checker, or the forcer.
- If the checker handles it but the forcer rejects it, decide: implement DAG-level support for table-record broadcasting in forcing, or explicitly reject at checker level with a clear compile-time error.
- Add an integration test that compiles, forces, and transpiles a `tab // rec` expression end-to-end.

---

## Priority 3: Binding cross-reference validation

**Driven by:** `spec.md §Bindings`

### 4a. Forward-reference detection

**Status:** The scope checker (`_check_scope` / `_walk_scope`) only detects duplicate bindings. Forward references (a binding referencing a later binding in the same block) cause the semantic evaluator to fail with "Undefined variable" at evaluation time with a poor error message.

**Spec requirement:** "Binding statements in a block can reference each other, but no recursion." The current implementation allows referencing prior bindings but has no explicit forward-reference or cycle detection.

**What to do:**
- Add a forward-reference validation pass in `check.py` that scans each block:
  - Collect all binding names in the block.
  - For each binding value expression, identify all referenced names.
  - If a referenced name belongs to a later binding in the same block, flag a forward-reference error.
- Add tests for:
  - Legal forward reference (should still be caught).
  - Simple non-recursive cross-reference.
  - Recursive cycle (should be caught).

---

## Priority 4: Spec-compliance hardening

### 4a. Non-saturating record bindings in DAG

**Driven by:** `dag.md §2.5` (record binding contract), `progress.md`

**Status:** The `dag.md` contract requires that any `record` binding in the final DAG must have explicit field structure and must not stand in for an unflattened direct step-call argument. Records that directly saturate task calls are already flattened during forcing. Non-saturating records may still appear in the DAG with incomplete field metadata.

**What to do:**
- Audit `force.py` to trace every `Record` value that reaches the DAG output.
- Verify each is either:
  - A non-saturating value that is not a direct step-call argument (OK, must have full field structure).
  - A saturating value that has been flattened into individual input bindings.
- If any record reaches the DAG without full field structure, add the missing fields or raise a compile-time error.
- Add a test that emits a DAG with a non-saturating record and verifies the record binding has complete `fields`.

### 4b. Universal mapped-port classification

**Driven by:** `dag.md §3.2` (MappedStep), `progress.md`

**Status:** `_mapped_step_bindings` in `force.py` populates `scatter` and `broadcast` arrays for every mapped step. The spec requires every input parameter to appear in exactly one of the two arrays.

**What to do:**
- Add a validation pass in `DAG.validate()` (or after forcing) that checks:
  - Every mapped step has both `scatter` and `broadcast` keys in `map`.
  - Every input name in the mapped step's `input_schema` appears in exactly one of the two arrays.
  - No input name appears in both.
- Add a test that forces a mapped step and verifies the classification is total and disjoint.

### 4c. Full optionality propagation

**Driven by:** `dag.md §2.3` (ParamSpec), `progress.md`

**Status:** `dag.Input` has `optional: bool`, serialized to DAG JSON. The `?` suffix from task annotations is detected in `force.py:_input()` via `typ.endswith('?')`. Coverage may be incomplete for inferred workflow inputs.

**What to do:**
- Audit every path in `force.py` that creates an `Input` object and verify optionality is detected:
  - `force_value` for `ir.Input` (line 102).
  - `_input` helper (line 877).
  - `_refine_input_metadata` (line 709).
- Add a test that compiles a task with `?`-marked inputs and verifies optionality propagates to the DAG output.

### 4d. Final-DAG merge freedom

**Driven by:** `dag.md §2.5` (merge contract), `progress.md`

**Status:** `_flatten_merge_value` in `force.py` flattens merge trees during DAG finalization. Tests exist for record merges and non-flattenable values. Edge cases with deeply nested merge trees at the output level may not be covered.

**What to do:**
- Add a test that constructs a deeply nested merge tree (e.g., `(a // b) // (c // d)`) and verifies the final DAG has no Merge bindings.
- Add a test that verifies the merging of records with overlapping and non-overlapping fields produces correct output field bindings.
- If any merge pattern cannot be flattened, ensure it produces a clear compile-time error (not a transpiler-side rejection).

---

## Dependency graph

```
P1a (record bindings) ─────── ✅ complete
P1b (nested fields) ───────── ✅ complete
P1c (merge outputs) ───────── ✅ complete

P2 (table updates) ────────── ✅ complete

P3 (binding validation) ───── ✅ complete

P4a (non-saturating records) ── ✅ complete
P4b (mapped-port validation) ── ✅ complete
P4c (optionality audit) ─────── ✅ complete
P4d (merge hardening) ───────── ✅ complete
```

---

## File change summary

| File | What changes |
|------|-------------|
| `ir/dag.py` | Added `_validate_mapped_ports` in `DAG.validate()` (lines 196-229) |
| `ir/force.py` | (P2) Input+Record rejection in Update handling (lines 133-150); (P1c) fix `_unwrap_non_record` to preserve `Field(source=StepCall)` values; (P4b fix) `_mapped_step_bindings` scatter/broadcast overlap fix |
| `ir/test_force.py` | P1c tests (merge record+field preservation, end-to-end no-Merge), P4d tests (deeply nested merge trees), P4c tests (optionality), P4a tests (non-saturating records), P2 tests (table update errors). Fixed pre-existing test bugs: `ir.Ref` → `ir.Name`, `ForceEnv` kwarg usage |
| `semantic/wf/check.py` | (P3) forward-reference detection in `_walk_scope`, new `_collect_name_refs` and `_walk_refs` methods |
| `semantic/wf/test_check.py` | (P3) forward-reference tests |
| `transpile/cwl/emit.py` | No changes needed (P1a, P1b, P1c already working) |
| `transpile/cwl/test_emit.py` | No changes needed (tests already exist for record bindings, nested fields) |
