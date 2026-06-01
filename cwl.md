# CWL Transpilation: Status and Implementation Plan

---

## 1. What is Implemented

The CWL transpiler lives at `swl/transpile/cwl/emit.py`. It reads DAG JSON and emits a packed CWL v1.0 document (`cwlVersion: v1.0`) with a `$graph` containing one `Workflow` node and one `CommandLineTool` per distinct task.

### 1.1 Task step → CommandLineTool

| DAG feature | CWL output | Test |
|---|---|---|
| `step.type == "task"` with embedded script | `CommandLineTool` + `InitialWorkDirRequirement` materializing `script.sh` | `test_transpile_function_workflow` |
| Input params (`file`/`str`/`int`/`float`) | Tool `inputs` with `type: File`/`string`/`int`/`float` | same |
| Output params with interpolation | `outputBinding.glob` as `$(inputs.var + '.ext')` expression | `test_output_glob_uses_cwl_expression` |
| CPU resource (`run.cpu`) | `ResourceRequirement.coresMin` | `test_transpile_function_workflow` |
| Memory resource (`run.memory`) | `ResourceRequirement.ramMin` (value in MB) | same |
| Docker image (`run.image`) | `DockerRequirement.dockerPull` | same |
| `run.time` | Silently dropped (CWL has no standard time limit) | none |

### 1.2 Workflow block

| DAG feature | CWL output | Test |
|---|---|---|
| Workflow inputs | `Workflow.inputs` with `type` and optional `doc` | `test_transpile_function_workflow` |
| Step calls with bindings | `Workflow.steps[]` with `in` array wiring inputs | same |
| Input binding (workflow input) | Source `#main/input_name` | same |
| Step output binding | Source `#main/step_id/output_name` | same |
| Literal binding | `default` value on step input port | `test_rejects_non_task_workflow_output` (indirect) |
| Field projection (step output) | Source `#main/step_id/field_name` | `_canonical_binding` handling |
| Field projection (input field) | Source `#main/input_name/field_name` | same |
| Workflow outputs | `Workflow.outputs` with `outputSource` and inferred type | `test_transpile_function_workflow` |
| Sub-workflow step | Embedded `Workflow` node in `$graph` | `test_imported_workflow_output_transpiles_as_workflow_step` |

### 1.3 Mapped steps (`map` → scatter)

| DAG feature | CWL output | Test |
|---|---|---|
| `MappedStep` without `group_by` | `scatter` on input ports + `scatterMethod: dotproduct` | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` |
| Mapped workflow | Scattered step referencing sub-workflow | `test_batch_mapped_workflow_emits_scattered_subworkflow` |
| Mapped lambda (generated sub-workflow) | Generated `Workflow` node + scattered call | `test_batch_mapped_simple_lambda_emits_generated_scattered_subworkflow` |
| Root partial `map` | Scattered sub-workflow from partial application | `test_root_partial_map_transpiles_as_scattered_subworkflow` |
| `map.source` via `input` | Each schema column gets a workflow input + scatter port | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` |
| `map.source` via `table` | Column bindings resolved to workflow input references | same (table branch) |

### 1.4 Rejected (clear error)

| Feature | Error message | Test |
|---|---|---|
| `map_by` (any `group_by`) | `CWL transpilation does not yet support map_by grouping` | `test_map_by_transpile_reports_explicit_grouping_not_supported` |
| Merge bindings in step inputs | `merge values are not supported` | `test_rejects_merged_task_input_binding` |
| Record bindings in step inputs | `record values are not supported` | (via `_step_input_error`) |
| Function-valued outputs | `function outputs are not supported` | (via `_workflow_output_error`) |
| Literal workflow outputs | `literal outputs are not supported` | `test_rejects_non_task_workflow_output` |
| `expr` interpolation in output paths | `Unsupported interpolation for CWL glob` | `test_rejects_output_expr_interpolation` |
| Merge/Record workflow outputs | `merge outputs are not supported` / `record outputs are not supported` | (via `_workflow_output_error`) |
| Merge bindings at deserialization | `Unsupported step binding during deserialization` | `test_rejects_merged_task_input_binding` |
| Nested `Field(Field(...))` | `field source ... is not supported` | (via `_step_input_error`) |

---

## 2. What is Not Implemented

### 2.1 `map_by` (grouped scatter)

**Status:** Explicitly rejected at `emit.py:195`.

**DAG signal:** `MappedStep` with `map.group_by` set.

**Reference:** `map_by_cwl.md` describes the full approach.

**What's needed:** Replace the rejection with a two-phase emission:
- **Phase 1 — Grouping:** An `ExpressionTool` with `InlineJavascriptRequirement` that reads the table columns as arrays, partitions rows by the key column, and writes per-group JSON files.
- **Phase 2 — Scatter:** A generated wrapper `CommandLineTool` (or `Workflow` for mapped workflows) that reads each group JSON and invokes the original function's script once per row, collecting outputs.

**Edge cases to handle:**
- The mapped function is a task vs. a workflow (different wrapper generation)
- The mapped function's inputs include array-typed params (arrays-of-arrays — reject initially)
- Large number of groups (thousands of JSON files)
- `expr` interpolation inside the wrapped function (must reject)

### 2.2 Record bindings

**Status:** Rejected at `emit.py:238`.

**DAG signal:** Binding with `{source: "record", fields: {...}}`.

**What's needed:** A record binding represents a constructed value (e.g., `{ bam: calls.bam, outbase: "merged" }`). CWL has no native record literal. Options:
- **Option A:** Flatten the record at compiler time: if the record directly saturates a task call, expand its fields into individual step inputs. This is already partially done by forcing (partial application produces flat inputs when the record saturates a task). For non-saturating records, emit an `ExpressionTool` that constructs the record from individual inputs.
- **Option B:** Emit an intermediate `ExpressionTool` step that takes each field as a separate input and produces the record as a JSON output. Downstream steps reference the JSON's field via `$(step.field)`.

### 2.3 Merge bindings

**Status:** Rejected. The DAG serializer cannot even represent merge bindings in step bindings (`dag.py:_binding_to_binding_dict` raises on Merge), and the CWL transpiler rejects them.

**DAG signal:** Binding with `{source: "merge", left: ..., right: ...}`.

**What's needed:** Merge bindings arise from `A | B | C` pipeline desugaring and from explicit `//` operators. In most cases the compiler's chain desugaring already produces flat bindings. Remaining case: explicit record update expressions that reach the DAG. Fix should be at the compiler level (in `force.py`), not in the transpiler:
- Add merge-flattening logic to forcing so that merge bindings are resolved to flat per-field bindings before DAG output.
- If a merge cannot be flattened (e.g., merging two non-overlapping records that feed a single port), reject at the DAG validation level with a clear message.

### 2.4 Nested field projections

**Status:** Rejected at `emit.py:242`.

**DAG signal:** `Field(Field(...), name)` — a field projection whose source is itself a field projection.

**What's needed:** Nested field projections (`x.y.z`) require the intermediate value to be resolved. CWL cannot chain field accesses natively. Options:
- **Option A:** Flatten at compiler time: if `x` is a step output and its type is known, resolve `x.y.z` to a direct step output reference.
- **Option B:** Emit an intermediate `ExpressionTool` that takes `x` as input and emits `x.y.z`.

### 2.5 `expr` interpolation in output paths

**Status:** Rejected at `emit.py:336-337`.

**DAG signal:** Interpolation part with `kind: "expr"` (e.g., `${outbase / 2}`).

**What's needed:** `expr` parts are arbitrary expressions, not simple variable references. CWL `$( ... )` expressions are JavaScript — they can represent arbitrary expressions if `InlineJavascriptRequirement` is declared. Options:
- **Prefix CWL requirement:** If the task or workflow has any `expr` interpolation, add `InlineJavascriptRequirement` and emit the expression as-is wrapped in `$( ... )`.
- **Keep rejecting:** Arbitrary SWL expressions may reference variables or operators that don't map to CWL's JS expression scope. Conservative approach: continue rejecting.

### 2.6 Time resource

**Status:** Silently dropped at `emit.py:175-183`.

**DAG signal:** `run.time.value` is present but `_resource_requirement` only reads `cpu` and `memory`.

**What's needed:** CWL v1.0 has no standard time limit in `ResourceRequirement`. Options:
- Emit as a hints entry: `{"class": "Hint", "timeLimit": <value>}`
- Emit in a non-standard extension namespace
- Continue dropping (current behavior)

### 2.7 Optional type serialization

**Status:** DAG does not serialize optionality. The `?` suffix in task annotations (e.g., `file?`) is used during semantic checking but is not present in the DAG output.

**What's needed:** If optionality is to be represented in CWL:
- Store optional flag in the DAG `InputSpec`
- CWL represents optional types as `["null", "File"]` (union with null)
- This depends on DAG schema changes first

### 2.8 Merge workflow outputs

**Status:** Rejected at `emit.py:255`.

**DAG signal:** `dag.outputs` containing `{source: "merge", ...}`.

**What's needed:** Same root cause as §2.3. Fix at compiler level: flatten merge bindings before DAG output so the transpiler never sees them.

### 2.9 Record workflow outputs

**Status:** Rejected at `emit.py:257`.

**DAG signal:** `dag.outputs` containing `{source: "record", ...}`.

**What's needed:** Records at workflow output level mean the workflow returns a constructed record, not a step output. Options:
- Decompose at compiler time: add synthetic steps to construct the record
- Support at transpiler time: emit `ExpressionTool` output construction

---

## 3. Gap Analysis Summary

| Feature | Severity | Root cause layer | Fix layer |
|---|---|---|---|
| `map_by` (grouped scatter) | **Blocking** for `map_by` workflows | Transpiler lacks two-phase emission | Transpiler: add ExpressionTool grouping + wrapper scatter (map_by_cwl.md) |
| Record bindings in step inputs | High | Transpiler rejects; compiler may not flatten all cases | Both: flatten in compiler, support remainder in transpiler |
| Merge bindings | High | Chain desugaring handles `\|`, explicit `//` may not flatten | Compiler: merge-flatten in force.py |
| Nested field projections | Medium | Compiler emits `Field(Field(...))` for multi-level access | Transpiler: ExpressionTool intermediary |
| `expr` interpolation | Medium | SWL arbitrary expressions don't map to CWL JS | Transpiler: optional InlineJavascriptRequirement + passthrough |
| Time resource | Low | CWL lacks standard time limit | Transpiler: emit as hints entry |
| Optional type serialization | Low | DAG schema missing optional flag | DAG schema: add optional field |
| Merge workflow outputs | Medium | Same as merge bindings | Compiler: merge-flatten |
| Record workflow outputs | Low | Rare in practice; compiler may emit for lambda results | Transpiler: ExpressionTool construction |
| `run` params missing in DAG | Low | Workflow steps' definition.run may be empty when inherited | Compiler: propagate run params to definition |

---

## 4. Implementation Plan

### Phase 1: `map_by` via ExpressionTool grouping + scattered wrapper

**Goal:** Replace `_validate_supported` rejection of `group_by` with full two-phase emission.

**Implementation steps:**

1. **Add `_validate_map_by_preconditions(step)`** — check that:
   - `map.source` is present and is either `input` or `table`
   - `map.group_by` names a column in `map.source.columns` or `step.input_schema`
   - The mapped function's inputs do not include array types (arrays-of-arrays not representable; reject with clear message)
   - The mapped function's output paths contain no `expr` interpolation

2. **Implement `_emit_grouping_expression_tool(step)`** — generate a CWL `ExpressionTool`:
   - Inputs: one per column in `step.input_schema` (typed as arrays matching the CWL array-of-items form), plus a `key_name: string` input
   - Requirements: `InlineJavascriptRequirement`
   - Outputs: `groups: File[]` with `outputBinding.glob: "groups/*.json"`
   - Expression (JavaScript): partition rows by key, write one JSON file per group using `fs.writeFileSync`
   - Structure of each group JSON:
     ```json
     { "key": "sample_A", "columns": { "fastq1": [...], "fastq2": [...] } }
     ```

3. **Implement `_emit_group_wrapper_tool(step)`** — generate a `CommandLineTool`:
   - Input: a single `File` (the group JSON); additional broadcast inputs from non-scatter bindings
   - Script: reads the group JSON, iterates over rows, calls the original tool's script for each row, collects outputs
   - Outputs: generated from `step.output_schema` (arrays — each group produces one output row)

4. **Implement `_emit_map_by_workflow_step(step, ...)`** that emits two CWL workflow steps:
   - Grouping step: calls the ExpressionTool
   - Processing step: calls the wrapper tool, `scatter` on the group file input

5. **Handle the mapped function being a sub-workflow** (not just a task):
   - The wrapper is a generated `Workflow` (using `SubworkflowFeatureRequirement`) that reads group JSON via inline `ExpressionTool`, calls the inner sub-workflow per row, and collects outputs.

6. **Update `_step_to_cwl` or call site:** When `step.map.group_by` is present, call `_emit_map_by_workflow_step` instead of the normal scatter path.

7. **Tests:**
   - `test_map_by_task_transpiles_to_grouping_and_scatter` — full DAG with simple `map_by` on a task
   - `test_map_by_workflow_transpiles` — `map_by` on a sub-workflow
   - `test_map_by_rejects_array_inputs` — verify rejection of array-typed inputs in mapped function
   - `test_map_by_rejects_expr_interpolation` — verify rejection of `expr` in wrapped function's outputs

**Complexity assessment:** ~250-350 lines of new Python code + ~100 lines of tests. The ExpressionTool JS and the wrapper bash script are the riskiest parts — they need careful quoting and edge-case handling.

### Phase 2: Merge flattening at compiler level

**Goal:** Eliminate merge bindings from DAG output so the transpiler never sees them.

**Implementation steps:**

1. **Identify merge sources:** Merge bindings in the DAG come from:
   - Pipeline chain desugaring (`A | B` → `x // a` bindings) — already flat in most cases
   - Explicit record update expressions (`r1 // r2`) in the SWL source
   - Partial application with `//` in the call site

2. **Add flattening pass in `force.py`** (after IR forcing, before DAG construction):
   - Walk all step bindings. When a binding is a `Merge`, recursively decompose it:
     - If both sides are `Input` or `Field(StepCall)`, flatten to individual per-field bindings
     - If one side cannot be decomposed (e.g., a `TableSource`), raise a compile-time error
   - Walk `dag.outputs`. Same treatment.

3. **Remove rejection in `_validate_supported`:** Once merges are flattened by the compiler, the transpiler never encounters them.

4. **Tests:**
   - Verify that existing pipeline tests (`function.swl`, `pipe.swl`) produce no merge bindings in DAG
   - Add test for explicit `//` in source that should be flattened
   - Add test for `//` that cannot be flattened (e.g., record-update on a table column) → expected compile error

### Phase 3: Record bindings via ExpressionTool

**Goal:** Support record bindings that reach the transpiler (non-saturating records).

**Implementation steps:**

1. **Determine if record can be saturated:** If the record binding is in a step's bindings for a named input, check whether the step's task expects a single record-typed input. If yes, the record is "saturating" and its fields should be flattened into individual inputs at compiler time.

2. **For non-saturating records at transpiler time:**
   - Generate an `ExpressionTool` that takes each record field as a separate input
   - The tool constructs the record as a JSON object and emits it as a `File`
   - Replace the record binding in downstream steps with a reference to the ExpressionTool's file output

3. **Tests:**
   - Create a DAG where a record does not saturate a task call
   - Verify that transpilation produces an intermediate ExpressionTool step

### Phase 4: Time resource hints

**Goal:** Preserve `run.time` in CWL output (currently dropped).

**Implementation steps:**

1. **Modify `_resource_requirement` (or create `_emit_runtime`):** After emitting `ResourceRequirement`, emit `time` as a hints entry:
   ```python
   hints = []
   if 'time' in run:
       hints.append({'class': 'TimeLimit', 'timeLimit': run['time'].get('value')})
   ```

2. **If hint approach is too non-standard:** Emit as a custom extension:
   ```python
   {'extension': 'https://swl-lang.org/timeLimit', 'value': run['time'].get('value')}
   ```

3. **Tests:**
   - Verify that function workflow CWL output includes time hint

### Phase 5: `expr` interpolation passthrough

**Goal:** Support `expr` interpolation parts by adding `InlineJavascriptRequirement` and passing the expression through.

**Implementation steps:**

1. **Modify `_interp_to_cwl_glob`:** Instead of raising on `expr` parts, emit the expression text as-is within `$( ... )`.

2. **Add `InlineJavascriptRequirement` to tool requirements** if any output uses `expr` interpolation.

3. **Tests:**
   - Create a `.swl` file with `expr` interpolation in output path
   - Verify transpilation succeeds and includes `InlineJavascriptRequirement`
   - Verify the emitted expression preserves the original expression text

### Phase 6: Nested field projections

**Goal:** Support `Field(Field(...), ...)` patterns.

**Implementation steps:**

1. **Add `_emit_field_resolution_expression_tool(binding)`:** For a nested field chain like `x.y.z`:
   - Create an ExpressionTool that takes the source binding as input
   - Resolve the field chain via JS expression: `$(inputs.source.y.z)`
   - Output the resolved value

2. **Wire the ExpressionTool output** as the input to the downstream step.

3. **Tests:**
   - Construct a DAG with `Field(Field(step_output, "inner"), "leaf")`
   - Verify transpilation produces an intermediate ExpressionTool

### Phase 7: Optional type support

**Goal:** Represent optional types in CWL as union types.

**Implementation steps:**

1. **Add optional flag to DAG `InputSpec`:** Requires DAG schema change (add optional `optional: bool` field).

2. **Modify `_cwl_type`:** When `optional=True`, emit `["null", <type>]` instead of just `<type>`.

3. **Tests:**
   - Transpile a DAG with optional input and verify CWL type is `["null", "File"]`

---

## 5. Effort Summary

| Phase | Feature | Layer | Estimated effort |
|---|---|---|---|
| 1 | `map_by` | Transpiler | ~350 lines emit.py, ~100 lines tests |
| 2 | Merge flattening | Compiler (force.py) | ~80 lines force.py, ~40 lines tests |
| 3 | Record bindings | Transpiler | ~150 lines emit.py, ~60 lines tests |
| 4 | Time resource hints | Transpiler | ~20 lines emit.py, ~20 lines tests |
| 5 | `expr` interpolation | Transpiler | ~30 lines emit.py, ~30 lines tests |
| 6 | Nested field projections | Transpiler | ~100 lines emit.py, ~50 lines tests |
| 7 | Optional types | DAG schema + transpiler | ~40 lines across dag.py + emit.py, ~30 lines tests |

**Total estimate:** ~800 lines Python + ~330 lines tests.

---

## 6. Dependencies and Ordering

```
Phase 2 (merge flattening) ─┐
                              ├──→ Phase 6 (field projections) can start
Phase 1 (map_by) ────────────┤     after Phase 2 (fewer edge cases)
                              │
Phase 3 (record bindings) ───┘

Phase 4 (time hints) ─── independent
Phase 5 (expr interp) ─── independent
Phase 7 (optional types) ─── depends on DAG schema change
```

Phase 2 should be done first because it cleans up the DAG and reduces the number of edge cases the transpiler must handle. Phase 1 is the highest user value (unblocks `map_by`). Phases 4-7 are incremental improvements.

---

## 7. Reference: Key Code Locations

| What | File | Lines |
|---|---|---|
| Transpiler entry point | `swl/transpile/cwl/emit.py` | 7-49 |
| Task → CommandLineTool | `_tool_to_cwl` | 52-89 |
| Step → CWL step (includes scatter) | `_step_to_cwl` | 111-141 |
| Workflow inputs/outputs | `_workflow_input_to_cwl`, `_workflow_output_to_cwl` | 92-108 |
| Binding validation (rejections) | `_validate_supported` | 193-211 |
| Canonical binding classification | `_canonical_binding` | 214-226 |
| CWL type mapping | `_cwl_type` | 306-323 |
| Interpolation → CWL glob | `_interp_to_cwl_glob` | 326-339 |
| `map_by` rejection | `_validate_supported` | 195-196 |
| `map_by` analysis | `map_by_cwl.md` | full file |
| DAG data model | `swl/ir/dag.py` | full file |
| DAG JSON format specification | `dag.md` | full file |
