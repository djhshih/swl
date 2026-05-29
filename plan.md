# Batch workflow / CWL generalization plan

## Goal
Make CWL emission work from a general canonical binding model instead of ad hoc shape checks, while keeping forcing as the single workflow definition builder.

## Problem observed
Nested/generated workflow helpers can currently reach CWL emission with raw internal binding shapes that the emitter does not understand, e.g. unresolved function outputs serialized as dict-like payloads. The emitter currently validates and emits by matching a small set of concrete Python classes, which is too narrow.

## Direction
1. Keep a single workflow definition builder in `python/swl/ir/force.py`.
2. Improve force/materialization so supported workflow bodies lower to concrete DAG values.
3. Make the CWL emitter more general by normalizing bindings to a canonical wireable subset before validation/emission.
4. Align the plan with the `tab` decision rather than any prior `[rec]` model.
5. Remove the separate `ArrayField` operator from parser/AST/IR/etc.; table column access should remain ordinary field access with type/schema-directed behavior.

## Canonical wireable binding subset
The emitter should normalize raw DAG bindings/outputs into one of:
- `input(name)`
- `step_output(step_id, output_name)`
- `array_step_output(step_id, output_name)`
- `input_field(input_name, field_name)`
- `literal(value)`

For batch workflows, `array_step_output(...)` should be understood as an array-valued/table-column-capable source arising from `tab` semantics, not from a separate array-of-record field operator.

Anything else should be either:
- normalized into one of those forms, or
- rejected before/at emission with a precise error.

## Emitter refactor
In `python/swl/transpile/cwl/emit.py`:
- add `_canonical_binding(value)`
- rewrite:
  - `_workflow_output_error()`
  - `_step_input_error()`
  - `_binding_source()`
  - `_infer_output_type()`
  to use `_canonical_binding()`
- keep workflow steps as first-class CWL `Workflow` runs in the `steps` block
- support input-field workflow outputs generically

## Force/materialization follow-up
In `python/swl/ir/force.py`:
- make `_materialize_workflow_dag()` guarantee that supported workflow outputs are concrete wireable values
- reject unresolved function outputs earlier with a clearer error if materialization still fails
- preserve enough type/schema information for `tab` field projection without introducing a dedicated `ArrayField` value kind unless forcing truly needs a serialization-only canonical form

## Parser / AST / IR alignment follow-up
Across parser, AST, semantic typing, IR, and DAG layers:
- remove `ArrayField` operator support from parser/AST/IR/etc.
- ensure `.` remains the single field-access operator
- make field access type-directed: record field on `rec`, column projection on `tab`
- reject missing table columns statically using schema information

## Why this is general
This avoids adding a special emitter case per function structure. Instead:
- forcing materializes workflow bodies
- emitter only needs to understand a small canonical set of wireable bindings
- nested workflows, generated helpers, and mapped workflows all go through the same emission model
- the same model works for `tab` without carrying a parallel array-of-record operator stack
