# Remaining implementation plan

## Goal

Drive `python/swl/` toward production readiness. Tracks concrete implementation work derived from:
- `dag.md` — the DAG JSON spec and cross-transpiler requirements
- `cwl.md` — CWL transpiler status and gaps
- `wdl.md` — WDL 1.1 transpiler implementation plan
- `nf.md` — Nextflow transpiler implementation plan
- `spec.md` — SWL language specification
- `issues.md` — prioritized issue tracker

---

## Priority 0: Merge flattening in forcing

### Why (dag.md §R2)

Merge bindings (`{source: "merge"}`) reach the DAG but are rejected by every transpiler — CWL raises `ValueError`, WDL raises `ValueError`, Nextflow has no clean mapping. Flattening must happen once in the compiler, not N times in N transpilers.

### Scope

- Walk step bindings after IR forcing, before DAG construction
- Recursively decompose `Merge` bindings into per-field flat bindings
- Walk `dag.outputs` — same treatment
- If a merge cannot be flattened (e.g., two non-overlapping record-valued step outputs into a single port), raise compile-time error

### Files

- `ir/force.py` — add flattening pass
- `ir/dag.py` — may need helper to detect merge-free DAG

### Tests

- Existing `test_rejects_merged_task_input_binding` must change from "rejected by transpiler" to "never appears in DAG"
- Add test for explicit `//` that should flatten
- Add test for `//` that cannot flatten → expected compile error

### After this step

- No `Merge` binding ever reaches a transpiler
- Each transpiler's merge rejection code can be removed
- `dag.md §2.5` merge flattening requirement is satisfied

---

## Priority 1: DAG metadata improvements

### 1a. Serialize `map.scatter` / `map.broadcast` in MappedStep

**Why (dag.md §R1):** Every transpiler independently re-derives which input ports scatter and which broadcast by inspecting `map.source.columns` and cross-referencing against `bindings`. This inference is fragile and duplicated.

**What to do:**
- Add `scatter` (array of port names) and `broadcast` (array of port names) to the `map` field of `MappedStep`
- Every input must appear in exactly one of the two lists
- Populate during forcing

**Files:** `ir/dag.py`, `ir/force.py`

### 1b. Serialize `optional` flag in InputSpec / ParamSpec

**Why (dag.md §R1):** Optionality (`file?`) is used during semantic checking but lost at DAG serialization. Transpilers need it for correct type declarations (CWL: `["null", "File"]`, WDL: `File?`, Nextflow: `optional`).

**What to do:**
- Add `optional: bool` field to `InputSpec` and `ParamSpec` (default: false)
- Preserve from task annotation parsing through to DAG serialization

**Files:** `ir/dag.py`, `syntax/task/parser.py` (ensure `?` flag is passed through)

### 1c. Serialize output types in top-level `outputs`

**Why (dag.md §R1):** Currently `outputs` is `{name: <Binding>}` with no type info. Every transpiler infers the type from the binding's source step — fragile when the source is an `Input`, `Literal`, or `Record`.

**What to do:**
- Change `outputs` from `{name: <Binding>}` to `{name: {type: <type>, desc: <desc>, value: <Binding>}}`
- Populate `type` from the inferred workflow output type during forcing

**Files:** `ir/dag.py`, `ir/force.py`, then each transpiler

---

## Priority 2: CWL transpiler — implement missing features

Driven by `cwl.md`. See that file for detailed implementation steps and test plans.

### 2a. `map_by` via ExpressionTool grouping + scattered wrapper (`cwl.md §4 Phase 1`)

**Why:** Currently `map_by` is rejected. This unblocks CWL users from using grouped scatter.

**What to do:**
- Implement `_emit_grouping_expression_tool(step)` — generates a CWL `ExpressionTool` with `InlineJavascriptRequirement` that partitions rows by key and writes per-group JSON files
- Implement `_emit_group_wrapper_tool(step)` — generates a `CommandLineTool` that reads a group JSON and calls the original tool for each row
- Implement `_emit_map_by_workflow_step(step)` — emits two CWL workflow steps (grouping + scatter)
- Handle both task and sub-workflow mapped functions

**Files:** `transpile/cwl/emit.py`, `transpile/cwl/test_emit.py`

### 2b. Record bindings via ExpressionTool (`cwl.md §4 Phase 3`)

**Why:** Record bindings are currently rejected. Non-saturating records need transpiler-side construction.

### 2c. Time resource hints, `expr` interpolation passthrough, nested fields (`cwl.md §4 Phases 4-6`)

Lower priority. Time can be emitted as hints. `expr` can pass through with `InlineJavascriptRequirement`. Nested fields can be resolved via intermediate ExpressionTool.

---

## Priority 5: Remaining spec-compliance gaps

### 5a. Table update semantics (`tab // rec`, `rec // tab`, `tab // tab`)

~~Currently fails with explicit error. Full implementation deferred (rarely used). Tighten error messages if needed.~~

**Status: DONE.** The semantic checker (`semantic/wf/check.py`) already handles `tab // rec`, `rec // tab`, and `tab // tab`. The force phase (`ir/force.py`) raises `ValueError` when a `StepCall` appears on either side of `//` because the forcing logic cannot flatten merges involving step call results. Error message tightened. Test added.

**Files:** `semantic/wf/check.py`, `ir/force.py`

### 5b. DAG circularity validation

~~Add validation pass in `DAG.validate()` or `force.py` that checks for cycles before emitting.~~

**Status: DONE.** `DAG.validate()` in `ir/dag.py` already implements DFS-based cycle detection and is called from `Forcer._finalize_dag()` in `ir/force.py`. Test for cycle detection added.

**Files:** `ir/dag.py` or `ir/force.py`

### 5c. `expr` interpolation validation

~~Currently `expr` parts are accepted by the compiler but rejected by CWL transpiler. The `dag.md §6.2` table clarifies per-target handling. Add per-transpiler validation.~~

**Status: DONE.** Per-transpiler validation already exists:
- **CWL:** `_interp_to_cwl_glob()` rejects `expr` parts with `ValueError` (tested in `test_rejects_output_expr_interpolation`)
- **WDL:** `_interp_to_wdl()` accepts `expr` parts as `~{text}` — WDL supports expressions in string interpolation
- **Nextflow:** `_interp_to_nf()` accepts `expr` parts as `${text}` — Nextflow supports shell-compatible expressions

---

## File change summary

| File | What changes |
|------|-------------|
| `ir/dag.py` | Add `scatter`/`broadcast` to map field; add `optional` to InputSpec/ParamSpec; change outputs format to include `type`/`desc`/`value` |
| `ir/force.py` | Add merge flattening pass; populate `map.scatter`/`map.broadcast`; populate output types; tightened table update error message (P5a) |
| `ir/test_force.py` | Added tests for table update error (P5a) |
| `ir/test_force_codec.py` | Added tests for DAG circularity validation (P5b) |
| `transpile/cwl/emit.py` | Implement `map_by` (Phase 1), record bindings (Phase 3), time hints (Phase 4), expr passthrough (Phase 5), nested fields (Phase 6); remove merge rejection |
| `transpile/cwl/test_emit.py` | Tests for all new features |
| `transpile/wdl/emit.py` | New file — WDL 1.1 transpiler |
| `transpile/wdl/test_emit.py` | New file — WDL transpiler tests |
| `transpile/nextflow/emit.py` | New file — Nextflow transpiler |
| `transpile/nextflow/test_emit.py` | New file — Nextflow transpiler tests |

---

## Dependency graph

```
P0 (merge flattening)
  └─→ P1a (scatter/broadcast ports) ──→ P2a (CWL map_by)
                                      ──→ P3 Phase 4 (WDL map_by)
                                      ──→ P4 Phase 4 (Nextflow map_by)
  └─→ P1b (optional serialization) ──→ P3 Phase 7 (WDL optional types)
  └─→ P1c (output types) ──────────→ P2, P3, P4 (all transpilers remove type inference)

P2 (CWL gaps) ─── independent, highest user value

P3 (WDL transpiler) ─── depends on P0, P1
P4 (Nextflow transpiler) ─── depends on P0, P1
P5 (spec gaps) ─── independent, incremental
```
