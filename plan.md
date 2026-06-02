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

**Status:** Record bindings in step inputs are accepted by CWL (`_step_input_error` returns `None`) but the code path is untested. Non-saturating records that reach the DAG may not produce correct CWL.

**What to do:**
- Add a DAG-level test that has a non-saturating record reaching CWL transpiler.
- If the record is directly passed to a step input, emit an intermediate `ExpressionTool` that constructs the record from individual input fields.
- Verify CWL output is valid.

### 1b. Resolve nested `Field(Field(...))` projections

**Status:** The CWL transpiler rejects `Field` with unsupported source types via `_step_input_error`. Simple field projections on step outputs (`Field(StepCall, name)`) work. Chains like `Field(Field(StepCall, "inner"), "leaf")` fail.

**What to do:**
- Flatten nested field chains at DAG construction time (in `force.py`) by resolving intermediate types to direct step-output references when possible.
- For cases where flattening is impossible (dynamic intermediate type), emit an `ExpressionTool` that takes the source and resolves the chain via JS expression.

### 1c. Merge workflow outputs

**Status:** `_workflow_output_error` rejects merge bindings at workflow output level. Flattening should happen in the compiler (`force.py`) but merge trees at output level may not be fully covered.

**What to do:**
- Audit `_final_outputs` / `_canonicalize_merges` / `_flatten_merge_value` to ensure all merge patterns at the output level are flattened.
- Add targeted tests for merge trees at output level that currently leak through.
- If a merge cannot be flattened, raise a compile-time error with a clear message.

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
P1a (record bindings) ─────── independent
P1b (nested fields) ───────── independent
P1c (merge outputs) ───────── dependent on P4d (merge hardening)

P2 (table updates) ────────── independent, low urgency

P3 (binding validation) ───── independent

P4a (non-saturating records) ── audit, likely low code change
P4b (mapped-port validation) ── new validation pass
P4c (optionality audit) ─────── audit, likely low code change
P4d (merge hardening) ───────── test coverage, may reveal edge cases
```

---

## File change summary

| File | What changes |
|------|-------------|
| `ir/dag.py` | (P4b) mapped-port validation in `DAG.validate()` |
| `ir/force.py` | (P1c) flatten remaining merge output edge cases; (P4a) ensure non-saturating records have full fields; (P4c) optionality audit |
| `ir/test_force.py` | Tests for P1c, P4a, P4c, P4d |
| `semantic/wf/check.py` | (P3) forward-reference validation pass |
| `semantic/wf/test_check.py` | Tests for P3 |
| `transpile/cwl/emit.py` | (P1a) record bindings via ExpressionTool; (P1b) nested field flattening or intermediate tool |
| `transpile/cwl/test_emit.py` | Tests for P1a, P1b, P1c |
| `compile.py` | (P4b) wire mapped-port validation into compilation pipeline |
