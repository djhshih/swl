# Batch workflow / CWL generalization plan

## Goal
Make CWL emission work from a general canonical binding model instead of ad hoc shape checks, while keeping forcing as the single workflow definition builder.

## Problem observed
Nested/generated workflow helpers can currently reach CWL emission with raw internal binding shapes that the emitter does not understand, e.g. unresolved function outputs serialized as dict-like payloads. The emitter currently validates and emits by matching a small set of concrete Python classes, which is too narrow.

## Direction
1. Keep a single workflow definition builder in `python/swl/ir/force.py`.
2. Improve force/materialization so supported workflow bodies lower to concrete DAG values.
3. Make the CWL emitter more general by normalizing bindings to a canonical wireable subset before validation/emission.

## Canonical wireable binding subset
The emitter should normalize raw DAG bindings/outputs into one of:
- `input(name)`
- `step_output(step_id, output_name)`
- `array_step_output(step_id, output_name)`
- `input_field(input_name, field_name)`
- `literal(value)`

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

## Why this is general
This avoids adding a special emitter case per function structure. Instead:
- forcing materializes workflow bodies
- emitter only needs to understand a small canonical set of wireable bindings
- nested workflows, generated helpers, and mapped workflows all go through the same emission model
