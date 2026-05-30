# Spec-compliance implementation plan

## Goal
Define the concrete code changes needed in `python/swl/` so the implementation matches the semantics in `spec.md` and the batch design in `new.md`.

This plan is intentionally code-focused. It does not restate the design at length; it identifies what to change, where to change it, and what behavior each change must enforce.

---

## Priority 1: make workflow typing explicitly represent `rec`, `tab`, and workflow function kinds

### Why
Right now workflow typing in `python/swl/semantic/wf/check.py` is mostly encoded through ad hoc evaluator values (`OpenRecord`, `ClosedRecord`, `TableValue`, `FunctionValue`, etc.) plus `TaskSignature` and `is_batch`. That is not enough to implement the spec cleanly.

The spec requires explicit distinction between:
- `rec`
- `tab`
- scalar arrays like `[file]`
- simple workflow: `rec -> rec`
- batch workflow: `tab -> rec` or `tab -> tab`

### Concrete changes

#### 1. Add workflow type definitions
Create a new module:
- `python/swl/semantic/wf/type.py`

Add workflow-level type objects, for example:
- scalar types: `file`, `str`, `int`, `float`
- array scalar types: `[file]`, `[str]`, `[int]`, `[float]`
- `RecordType(fields, open=True/False)`
- `TableType(columns)`
- `FunctionType(input_type, output_type)`
- helpers like `is_record_type`, `is_table_type`, `is_function_type`

This layer must be independent from `TaskSignature`, though task signatures should be convertible into workflow `FunctionType(rec, rec)`.

#### 2. Refactor workflow semantic checking to use workflow types
Modify:
- `python/swl/semantic/wf/check.py`

Changes required:
- replace `batch=True` / `is_batch=True` as the main classifier with workflow-type inference
- infer whether lambda parameter is `rec` or `tab` from how it is used
- infer field requirements on records and columns on tables
- represent `map f` as taking `rec -> rec` and returning `tab -> tab`
- reject `map` when callee is batch-typed
- reject missing table columns statically
- classify workflow root by inferred function type, not by `_uses_map()` heuristics alone

Do not remove all evaluator-style helper values immediately if not necessary, but they must become an implementation detail behind explicit type inference.

#### 3. Expose workflow type info in semantic results
Extend `WorkflowCheck` in:
- `python/swl/semantic/wf/check.py`

Add fields such as:
- `workflow_type`
- `root_input_type`
- `root_output_type`

Keep `signature` if useful for compatibility, but it should become a projection from the workflow type, not the canonical source of truth.

### Required behavior after this step
- simple workflows are explicitly classified as `rec -> rec`
- batch workflows are explicitly classified as `tab -> rec` or `tab -> tab`
- `map` type-checking is based on function types
- missing `tab` field access is rejected at semantic-check time

---

## Priority 2: bring builtins in line with the spec surface

### Why
The spec defines both `map` and `map_by`. Code implements only `map`.

### Concrete changes

#### 1. Implement explicit unsupported handling for `map_by`
Modify:
- `python/swl/semantic/wf/check.py`
- `python/swl/ir/lower.py`

Add builtin recognition for `map_by`, but initially reject it with a precise error such as:
- `map_by is not implemented`

This is required for honest spec compliance handling: unsupported spec features must fail clearly rather than behaving as unknown identifiers or malformed applications.

#### 2. Keep builtin lowering structured
In:
- `python/swl/ir/lower.py`

Refactor builtin detection so both `import`, `map`, and `map_by` use a single builtin-dispatch path.

### Required behavior after this step
- `map` continues to work
- `map_by` is recognized and rejected explicitly until implemented

---

## Priority 3: make IR semantics match the spec cleanly

### Why
The spec and new design require:
- `map` as an explicit semantic construct
- `.` as the only field operator
- no `ArrayField`

The current IR is close, but type/schema information is still too implicit.

### Concrete changes

#### 1. Keep `ir.Map`, keep `ir.Field`, do not add `ArrayField`
Modify if needed:
- `python/swl/ir/node.py`

Add optional type/schema fields to IR nodes where useful, especially:
- `Field` should be able to carry whether its source is record-like or table-like if semantic typing already knows that
- `Map` should be able to carry or reference the mapped callable’s input/output schema if needed later in forcing

Do not introduce a new field-access IR node.

#### 2. Lower builtin applications in a typed way
Modify:
- `python/swl/ir/lower.py`

When lowering:
- `map f xs` should lower to `ir.Map(function, arg)`
- ordinary `.` stays `ir.Field(...)`
- builtin info from semantic checking should be preserved enough that force/emitter do not have to rediscover all semantics heuristically

### Required behavior after this step
- IR continues to represent mapped execution explicitly
- no parallel field-access operator hierarchy exists

---

## Priority 4: separate canonical SWL DAG construction from CWL-specific flattening

### Why
This is the largest code/spec mismatch.

The spec says the workflow-language batch value is `tab`.
Current forcing partially flattens mapped `tab` input into per-column ports too early, inside `force.py`. That is backend lowering logic, not language semantics.

### Concrete changes

#### 1. Make `force.py` produce a canonical logical DAG
Modify:
- `python/swl/ir/force.py`
- `python/swl/ir/dag.py`

Change the canonical forced representation so it preserves:
- normal step call
- mapped step call
- mapped source as logical `tab` source
- ordinary field projection from mapped outputs
- task/workflow callee kind
- input/output schema for mapped callees

Specifically, revise:
- `Forcer._mapped_step_bindings()`
- `Forcer._tab_column_input_spec()`
- `Forcer._prune_unused_inputs()`

Current behavior to remove from canonical forcing:
- synthesizing flattened workflow inputs like `fastq1`, `fastq2`, etc. as the primary representation of the mapped source

Replace with:
- logical map source remains `xs` (or another table source)
- mapped step stores schema/signature metadata needed for later backend normalization

#### 2. Extend DAG classes to carry schema explicitly
Modify:
- `python/swl/ir/dag.py`

Extend `MappedStep` with explicit metadata such as:
- logical source kind
- input schema
- output schema
- callee kind already exists as `type`

Also make serialization/deserialization preserve this metadata.

### Required behavior after this step
- `eval_force.py` on batch workflows emits a symbolic DAG that still speaks in logical `tab` terms
- mapped execution remains symbolic
- canonical DAG does not depend on backend scatter-port flattening

---

## Priority 5: move scatter-port expansion into CWL normalization

### Why
CWL needs flattened scatterable ports; SWL semantics do not. That conversion must happen in the CWL layer.

### Concrete changes

#### 1. Add a CWL normalization step for mapped/table flow
Modify:
- `python/swl/transpile/cwl/emit.py`

Add functions that normalize canonical logical DAG bindings into CWL-ready sources, for example:
- logical mapped input source -> scattered input ports
- logical table column projection -> array-valued step output source

The current `emit.py` logic around:
- `_step_to_cwl()`
- `_canonical_binding()`
- `_binding_source()`
- `_infer_output_type()`

should be refactored so it works from canonical logical DAG semantics, not from forcing-side flattening.

#### 2. Generate scatter ports from mapped callee input schema
In `emit.py`, when emitting a mapped step:
- inspect mapped step’s input schema
- create one scatter input per mapped input field
- wire them from the normalized workflow-level representation
- use `scatterMethod: dotproduct`

#### 3. Keep reduction wiring from projected mapped outputs
Downstream inputs like:
- `merge { bam: calls.bam }`

must still emit as array-valued CWL sources.

### Required behavior after this step
- batch workflows transpile to runnable CWL using native scatter
- flattening occurs only in the CWL layer
- canonical DAG remains SWL-oriented

---

## Priority 6: implement explicit validation/rejection for unsupported table semantics

### Why
The spec defines more `tab` semantics than the code currently supports, especially around validation and `//`.

### Concrete changes

#### 1. Reject unsupported table-update semantics explicitly
Modify:
- `python/swl/semantic/wf/check.py`
- `python/swl/ir/force.py`

For now, if workflow expressions use:
- `tab // tab`
- `tab // rec`
- `rec // tab`

outside the limited supported cases, raise explicit errors such as:
- `table update semantics are not implemented`

Do not silently treat these as ordinary record merges if that would violate the spec.

#### 2. Prepare input validation hooks for concrete `tab` values
Modify or add:
- possibly a new runtime/pre-run validation module under `python/swl/semantic/wf/` or `python/swl/ir/`

If full external input validation is not implemented now, at minimum define one validation entry point that will eventually check:
- all table columns are arrays
- all arrays have equal length

If not wired yet, keep it stubbed but explicit.

### Required behavior after this step
- unsupported table semantics fail clearly
- there is no false impression that full spec `tab` update/validation is already implemented

---

## Priority 7: align compile/CLI surfaces with the new canonical artifact

### Why
The command-line entry points should reflect the canonical model after the refactor.

### Concrete changes

Review and adjust if needed:
- `python/swl/compile.py`
- `python/swl/eval_force.py`
- `python/swl/eval_ir.py`

Required changes:
- ensure `compile_workflow()` writes the canonical symbolic DAG
- ensure debug output from `eval_ir.py` and `eval_force.py` clearly shows mapped steps and logical sources
- ensure no CLI tool assumes all batch structure has already been flattened to backend ports

---

## Priority 8: update tests to enforce spec-compliant behavior

### Why
The current tests cover many good cases, but several of them encode the current premature-flattening behavior. Those tests should move to the CWL layer, while force/DAG tests should assert canonical SWL semantics.

### Concrete changes

#### 1. Semantic tests
Modify/add in:
- `python/swl/semantic/wf/test_check.py`

Add tests for:
- explicit workflow root type classification: `rec -> rec`, `tab -> rec`, `tab -> tab`
- missing `tab` column rejection
- `map_by` explicit unsupported error
- unsupported table-update forms rejected explicitly

#### 2. IR/force tests
Modify/add in:
- `python/swl/ir/test_force.py`
- `python/swl/ir/test_force_codec.py`

Change expectations so force tests assert:
- mapped steps remain symbolic
- canonical DAG preserves logical table source
- canonical DAG does not flatten `xs` into workflow-level column inputs as its primary meaning
- no legacy `ArrayField` encoding appears

#### 3. CWL tests
Modify/add in:
- `python/swl/transpile/cwl/test_emit.py`

Move flattening-specific assertions here:
- scattered per-column input ports are created during CWL emission
- `scatterMethod: dotproduct` is emitted
- mapped workflow/tool steps still run correctly
- reduction tasks consume array-valued mapped outputs

### Required behavior after this step
- tests reflect the proper layering:
  - semantic typing
  - canonical DAG
  - backend CWL lowering

---

## Priority 9: code-specific file checklist

### `python/swl/semantic/wf/type.py`
- new file
- add workflow type definitions and helpers

### `python/swl/semantic/wf/check.py`
- replace heuristic batch classification with explicit workflow typing
- add `map_by` recognition and explicit rejection
- add explicit unsupported errors for unimplemented table-update semantics
- expose root workflow type information

### `python/swl/semantic/task/type.py`
- leave task annotation typing intact
- add conversion helpers if needed from task signature to workflow function type

### `python/swl/ir/node.py`
- keep `Map`
- keep `Field`
- optionally attach inferred type/schema metadata
- do not add `ArrayField`

### `python/swl/ir/lower.py`
- route builtins through a single builtin-dispatch path
- preserve semantic type/schema info into IR where needed
- recognize `map_by` and lower to explicit unsupported handling or placeholder builtin node

### `python/swl/ir/dag.py`
- extend `MappedStep` to preserve logical source + schema metadata
- update `to_dict()` / `from_dict()` accordingly

### `python/swl/ir/force.py`
- stop performing CWL-style flattening in canonical DAG construction
- preserve logical `tab` map sources and mapped schemas
- explicitly reject unsupported table-update semantics

### `python/swl/transpile/cwl/emit.py`
- add normalization from canonical logical `tab` DAG to scatterable CWL ports
- refactor `_canonical_binding`, `_binding_source`, `_infer_output_type`, `_step_to_cwl`
- keep mapped output columns as array-valued sources in CWL

### `python/swl/compile.py`
- ensure canonical DAG is what gets written

### `python/swl/eval_ir.py`, `python/swl/eval_force.py`
- keep debug printing aligned with the new logical representation

### `python/swl/semantic/wf/test_check.py`
- add root-type and unsupported-feature tests

### `python/swl/ir/test_force.py`, `python/swl/ir/test_force_codec.py`
- update canonical DAG assertions

### `python/swl/transpile/cwl/test_emit.py`
- keep scatter flattening assertions here, not in force tests

---

## Recommended implementation order

1. add workflow type module
2. refactor semantic checker to infer explicit workflow/root types
3. add explicit `map_by` unsupported handling
4. refactor force/DAG so canonical batch representation stays logical
5. refactor CWL emitter to normalize logical `tab` flow into scatter ports
6. add explicit rejection for unsupported table-update semantics
7. update CLI/debug surfaces
8. update tests to match the new layering

---

## Definition of done

The codebase is compliant enough with the spec/new design when all of the following are true:
- workflow typing explicitly distinguishes `rec`, `tab`, and workflow function kinds
- `map` is typed as `(rec -> rec) -> tab -> tab`
- `.` remains the only field operator for both record fields and table columns
- batch workflows are represented canonically as `tab` workflows, not defined by flattened backend ports
- mapped execution remains symbolic in the compiled DAG
- CWL flattening occurs during CWL lowering, not during canonical DAG construction
- `map_by` is either implemented or explicitly rejected as unsupported
- unsupported table-update semantics fail clearly
- tests enforce the above separation and behavior