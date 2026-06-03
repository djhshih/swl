# Refactoring plan for `python/`

## Goals

Reduce redundancy, dead code, and incidental complexity while preserving the compiler pipeline described in `spec.md`, `dag.md`, and `architecture.md`:

- syntax/task parsing
- semantic checking
- IR lowering
- IR forcing
- normalized DAG emission
- transpilation from normalized DAG

Primary focus areas:

1. remove duplicated logic and stale compatibility shims
2. tighten the final DAG contract so transpilers stop compensating for compiler variability
3. split oversized modules by responsibility
4. delete debug artifacts and dead code
5. improve testability and maintenance boundaries

---

## Executive summary

The codebase already has a clear conceptual pipeline, but implementation-wise it is concentrated in a few large, mixed-responsibility files:

- `python/swl/semantic/wf/check.py` (1155 lines)
- `python/swl/ir/force.py` (1136 lines)
- `python/swl/transpile/cwl/emit.py` (742 lines)
- to a lesser extent `python/swl/ir/lower.py` (442 lines) and `python/swl/ir/dag.py` (474 lines)

The biggest maintainability problem is not just file size; it is **contract duplication**:

- type normalization appears in multiple places
- mapped-step semantics are partly in the forcer and partly re-inferred in transpilers
- bindings have **three** serialized shapes (`_binding_to_dict` in dag.py, `_binding_to_binding_dict` in dag.py, `_binding_to_public_dict` in force.py)
- workflow/task loading and validation logic is duplicated between checker and forcer
- transpilers still contain fallback logic for DAG forms that `dag.md` says should already be normalized away

The highest-value refactor is therefore:

> Move more normalization into the compiler core, make the DAG representation stricter and singular, and simplify transpilers into mostly direct renderers.

---

## Codebase findings

### 1. Oversized modules with mixed responsibilities

#### `python/swl/semantic/wf/check.py`
This file handles:

- import loading
- file reading
- bash variable validation
- scope validation
- chain type checking
- semantic evaluation
- input inference
- workflow signature building
- table/map/map_by semantics
- imported workflow execution modeling

This is too much in one module. It mixes:

- front-end validation
- abstract interpretation
- import infrastructure
- reusable utility logic

#### `python/swl/ir/force.py`
This file handles:

- forcing IR values
- environment management
- step emission
- workflow materialization
- task parsing/type extraction for embedded tool defs
- merge flattening
- output type inference
- DAG finalization
- binding serialization helpers

This is the compiler's most overloaded file.

#### `python/swl/transpile/cwl/emit.py`
This file handles:

- DAG validation (overlaps with `dag.py:DAG.validate()`)
- CWL rendering
- nested workflow packaging
- ad hoc record tool synthesis
- map and map_by lowering
- output type inference fallback
- binding canonicalization helpers (reverse-engineering DAG serialization)

It is doing compiler cleanup that should largely happen before transpilation.

---

### 2. Duplicated normalization logic

Examples:

- array-type wrapping appears in both `ir/force.py` and `transpile/cwl/emit.py` (`_as_array_type`)
- output type inference exists in `ir/force.py`, `transpile/cwl/emit.py`, and `transpile/wdl/emit.py`
- binding classification/canonicalization exists separately in `ir/dag.py` serialization code and `transpile/cwl/emit.py`
- workflow/task file parsing and bash-variable validation is performed in both `semantic/wf/check.py` and `ir/force.py`
- map input classification is produced in forcing, but CWL still has stale logic around `map.get('ports')`
- field-chain traversal duplicated: `dag.py:_nested_field_to_dict` and `cwl/emit.py:_field_chain` walk Field chains identically

This redundancy makes the code fragile when the DAG contract evolves.

---

### 3. Stale or dead code paths

Concrete examples found:

- `python/swl/compile.py`: `UserError` is defined but never raised by the current pipeline
- `python/swl/semantic/task/type.py`: `TypeChecker.infer_workflow_input()` appears unused
- `python/swl/ir/lower.py`: `_record_field_names()` appears unused
- `python/swl/syntax/wf/parser.py`: `_until()` and `_find()` appear unused
- `python/swl/ir/force.py`: duplicate `if isinstance(value, StepCall): return set(value.outputs)` in `_available_inputs` (lines 302-303 and 304-305)
- `python/swl/transpile/cwl/emit.py`: unreachable `return` on line 443 in `_cwl_type` after line 442 returns
- debug/evaluation scripts in package root not part of production pipeline:
  - `python/swl/eval.py`
  - `python/swl/eval_ir.py`
  - `python/swl/eval_force.py`
  - `python/swl/eval_task.py`
  - `python/swl/eval_task_semantic.py`
  - `python/swl/eval_wf_semantic.py`
- checked-in caches/artifacts:
  - `python/.pytest_cache/`
  - `python/swl/**/__pycache__/`
  - compiled `.pyc` files under `python/swl/**/__pycache__/`

These are immediate cleanup candidates.

---

### 4. Inconsistent DAG/binding surface

`dag.md` specifies a normalized DAG contract, but the implementation has **three** serialization formats for bindings:

| Function | Location | Step-output encoding | Used for |
|----------|----------|---------------------|----------|
| `_binding_to_dict` | `ir/dag.py:287` | `{'step': id, 'output': name}` | outputs, nested bindings |
| `_binding_to_binding_dict` | `ir/dag.py:353` | `{'source': id, 'output': name}` | step bindings in `to_dict()` |
| `_binding_to_public_dict` | `ir/force.py:1121` | `{'source': 'step', 'step': id, 'output': name}` | map source metadata in forcing |

Additionally:

- step bindings use a different mini-schema than general bindings (`_binding_to_binding_dict` omits `source` key for trivial Input bindings, returning `{}`)
- `transpile/cwl/emit.py` has to reverse-engineer "canonical bindings" via `_canonical_binding()` instead of receiving one normalized binding API
- `DAG.from_dict()` contains compatibility logic for heterogeneous encodings
- `_nested_field_to_dict` (dag.py:369) and `_field_chain` (cwl/emit.py:332) duplicate field-chain flattening
- transpilers still reject `Merge` and `Record` forms that `dag.md` says should be flattened before emission

This is a major source of bloat.

---

### 5. Dynamic pseudo-types and ad hoc objects

Example:

- `transpile/cwl/emit.py` creates ad hoc input specs with `type('InputSpec', (), {...})()` (lines 57, 62)

This is a maintainability smell. The project already uses dataclasses elsewhere; fake runtime types should be replaced with real structs/helpers.

---

### 6. Repeated stringly-typed type conversion

There is no single source of truth for:

- SWL scalar/array type names
- optional type handling
- conversion to CWL/WDL/NF types
- conversion from semantic types to DAG types

Currently this is spread across:

- `semantic/task/type.py`
- `semantic/wf/type.py`
- `ir/force.py` (`_normalize_swl_type`, `_as_array_type`, `_is_optional_type`)
- `transpile/cwl/emit.py` (`_cwl_type`, `_as_array_type`)
- `transpile/wdl/emit.py` (`_wdl_type`)
- `transpile/nf/emit.py` (`_input_qualifier`)

This invites drift.

---

### 7. Load/cache duplication

Import/file loading and definition caching are spread across checker, lowerer, and forcer:

- `semantic/wf/check.py:_load_imports` reads and parses source files
- `ir/force.py:_tool_definition` reparses `.sh` files and recomputes signatures
- `ir/lower.py:_cached_workflow_body` maintains its own cache of lowered workflows
- each stage knows too much about path resolution and relative imports

---

## Refactoring plan

## Phase 0 — Repository hygiene and deletion of obvious dead weight

### Actions

1. Remove generated artifacts from source tree:
   - `python/.pytest_cache/`
   - `python/swl/**/__pycache__/`
   - `*.pyc`

2. Add or tighten ignore rules for these artifacts.

3. Either remove or move debug scripts out of the package root:
   - `python/swl/eval.py`
   - `python/swl/eval_ir.py`
   - `python/swl/eval_force.py`
   - `python/swl/eval_task.py`
   - `python/swl/eval_task_semantic.py`
   - `python/swl/eval_wf_semantic.py`

   Recommended destination:
   - `python/tools/` or `python/devtools/`
   - or delete if they are no longer used

4. Remove clearly unused functions/classes:
   - `compile.py:UserError` — also audit whether the pipeline should raise it or remove the catch
   - `semantic/task/type.py:TypeChecker.infer_workflow_input`
   - `ir/lower.py:_record_field_names`
   - `syntax/wf/parser.py:_until`
   - `syntax/wf/parser.py:_find`

5. Fix dead/duplicate branches:
   - duplicate `StepCall` branch in `ir/force.py:_available_inputs` (lines 302-303 and 304-305)
   - unreachable second return in `transpile/cwl/emit.py:_cwl_type` after line 442

### Outcome

Low-risk cleanup that shrinks noise immediately and makes the next phases easier.

---

## Phase 1 — Define strict internal contracts first

This phase is the most important. It aligns with `dag.md` §R1 (explicit over inferred) and §R2 (flatten before emitting).

### 1.1 Create a single normalized binding model

Introduce a dedicated binding module, e.g.:

- `python/swl/ir/binding.py`
  - dataclasses or tagged classes for final DAG bindings
  - serialization/deserialization helpers
  - validation helpers

Consolidate binding-related logic out of scattered locations:

- from `ir/dag.py` (consolidate `_binding_to_dict`, `_binding_to_binding_dict`, `_nested_field_to_dict`)
- from `ir/force.py` (eliminate `_binding_to_public_dict`)
- from `transpile/cwl/emit.py` (eliminate `_canonical_binding`, `_field_chain`)

Target API:

- one in-memory binding model
- one JSON encoding for that model
- one deserializer
- one validator

### 1.2 Eliminate legacy/alternate encodings

Refactor `DAG.to_dict()` / `DAG.from_dict()` to emit and accept a single canonical format.

In particular:

- stop having separate serialization conventions for step bindings vs outputs
- stop using mixed forms like:
  - bare `{'step': ..., 'output': ...}`
  - `{'source': step_id, 'output': ...}`
  - `{'kind': 'field', ...}`
  - `{'source': 'step', ...}`
  - empty dict `{}` as a degenerate binding shape
- choose one encoding per binding kind and use it everywhere

Recommended JSON forms (per `dag.md`):

- input: `{ "source": "input", "name": ... }`
- step output: `{ "source": "step_output", "step": ..., "output": ... }`
- literal: `{ "source": "literal", "value": ... }`
- field: `{ "source": "field", "field": ..., "value": ... }`
- record: `{ "source": "record", "fields": ... }`
- table: `{ "source": "table", ... }`

Then deprecate support for alternate shapes.

### 1.3 Move DAG invariant checks into compiler core

Strengthen finalization so the final DAG always satisfies `dag.md` §R2 and the Binding contracts.

Compiler should guarantee before emission:

- no `Merge` in final DAG (per dag.md Binding contract: "Merge bindings must not appear in the final DAG")
- direct-call `Record` saturation completed (per dag.md Record contract: "If a record literal directly saturates a known task or workflow interface, the compiler must flatten it")
- normalized mapped-step metadata complete (`scatter`, `broadcast`) — every input in exactly one category
- explicit output types present (per dag.md OutputSpec contract: "Every emitted workflow output must have an explicit scalar or gathered-array type")
- explicit optionality preserved (per dag.md §R1: "Optional vs required parameters" currently lost at serialization)
- no transpiler-specific fallback assumptions needed

Once these are guaranteed, remove redundant validation from transpilers:
- `transpile/cwl/emit.py:_validate_supported` overlaps with `dag.py:DAG.validate()`
- `transpile/cwl/emit.py:_step_input_error` rejects Merge/Record/ForcedFunction — should be unreachable
- `transpile/cwl/emit.py:_workflow_output_error` rejects Merge/ForcedFunction — should be unreachable

### Outcome

This creates a clear compiler/transpiler boundary and removes the root cause of most transpiler complexity.

---

## Phase 2 — Split `ir/force.py` by responsibility

Current `ir/force.py` should be decomposed into smaller units. The exact split should serve the forcing stage's natural sub-tasks: context management, value evaluation, step emission, DAG finalization, tool definition loading, and merge handling.

### Proposed split

#### `python/swl/ir/forcing/context.py`
- `ForceEnv`
- shared force-time context/caches

#### `python/swl/ir/forcing/evaluator.py`
- `force_value`
- `_force_ref`
- `_force_apply`
- `_apply`
- lambda forcing

#### `python/swl/ir/forcing/emit.py`
- `_emit_task_call`
- `_emit_workflow_call`
- `_emit_mapped_step`
- `_mapped_step_bindings`
- step dependency helpers

#### `python/swl/ir/forcing/finalize.py`
- `_refine_input_metadata`
- `_flatten_step_bindings`
- `_flatten_outputs`
- `_build_output_specs`
- `_prune_unused_inputs`
- `_assert_wireable_output`

#### `python/swl/ir/forcing/tooldefs.py`
- `_tool_definition`
- `_workflow_definition`
- `_materialize_workflow_dag`
- task/workflow definition caching

#### `python/swl/ir/forcing/merge.py`
- merge canonicalization helpers
- flatten helpers
- output normalization helpers

### Additional cleanup in this phase

1. Stop reparsing task files in forcing if the checker already loaded them.
   - Today `_tool_definition()` reparses source via `TaskParser()` and recomputes signatures.
   - Instead carry imported task/workflow definitions forward from checker/lowerer.

2. Replace ad hoc helper functions with focused utility modules.

3. Add a small public façade:
   - `force_file()`
   - `Forcer.force()`

Everything else should become private to submodules.

### Outcome

Forcing becomes understandable and testable in pieces.

---

## Phase 3 — Split `semantic/wf/check.py` into front-end services

This file is too large and conceptually dense.

### Proposed split

#### `python/swl/semantic/wf/imports.py`
- `Import`
- file loading
- recursive workflow/task import handling
- circular import tracking

#### `python/swl/semantic/wf/scope.py`
- `_check_scope`
- reference collection
- duplicate/forward-reference logic

#### `python/swl/semantic/wf/infer.py`
- symbolic value classes:
  - `OpenRecord`
  - `ClosedRecord`
  - `FunctionValue`
  - `ClosureValue`
  - `ComputationValue`
  - `TableValue`
  - `UnknownValue`
  - `TypedValue`
- expression evaluation for inference
- map/map_by symbolic semantics

#### `python/swl/semantic/wf/signature.py`
- `_build_workflow_signature`
- helpers converting inferred values to `TaskSignature` and `wf_type`

#### `python/swl/semantic/wf/bashvars.py`
- `_validate_bash_variables`
- interpolation variable extraction

#### `python/swl/semantic/wf/check.py`
- keep only orchestration API:
  - `Checker.load()`
  - `Checker.load_content()`
  - assembly of results

### Additional cleanup

1. Convert semantic value classes to dataclasses where helpful.
2. Give them clearer names if needed (`SemanticOpenRecord`, etc.) to distinguish them from DAG/IR records.
3. Stop exposing too many internals via underscored methods that are reused elsewhere.

### Outcome

The semantic layer becomes a set of smaller, named subsystems rather than one monolith.

---

## Phase 4 — Centralize type and optionality handling

Create a single type utility module, e.g.:

- `python/swl/types.py`
- or `python/swl/semantic/types.py`

### Responsibilities

- normalize SWL type strings
- detect optional vs array
- strip/add optional markers
- strip/add array wrappers
- convert semantic task types to DAG types
- target-type conversion helpers per backend

### Example functions

- `normalize_swl_type(s: str | None) -> str | None`
- `is_optional_type(s)`
- `is_array_type(s)`
- `array_item_type(s)`
- `to_array_type(s)`
- `base_scalar_type(s)`
- `to_cwl_type(s)`
- `to_wdl_type(s)`
- `to_nf_qualifier(s)`

### Code to migrate

From:

- `semantic/task/type.py`
- `semantic/wf/type.py`
- `ir/force.py:_normalize_swl_type`, `_as_array_type`, `_is_optional_type`
- `transpile/cwl/emit.py:_cwl_type`, `_as_array_type`
- `transpile/wdl/emit.py:_wdl_type`
- `transpile/nf/emit.py:_input_qualifier`

### Outcome

Less stringly-typed drift and easier spec evolution.

---

## Phase 5 — Simplify transpilers after DAG cleanup

Once the DAG is stricter, simplify all transpilers to consume it directly.

### 5.1 CWL transpiler

Refactor `transpile/cwl/emit.py` into:

- `validate.py` for target-only checks (remove checks already guaranteed by compiler)
- `bindings.py` for CWL binding expression rendering (replace `_canonical_binding`)
- `tools.py` for `CommandLineTool`/subworkflow emission
- `workflow.py` for top-level workflow emission
- `map.py` for map/map_by-specific CWL logic

### Specific simplifications

- remove `_canonical_binding()` — normalized bindings already provide one shape
- remove `_field_chain()` — field chains handled by normalized binding model
- remove fallback output type inference (`_infer_output_type`) where `OutputSpec.type` is guaranteed
- remove stale `ports` logic and rely on `map.scatter` / `map.broadcast`
- replace fake dynamic `InputSpec` objects with real helper dataclasses or plain dict helpers
- remove `_validate_supported` checks that duplicate DAG-level invariants (Merge rejection, unsupported binding types)
- remove `_step_input_error` / `_workflow_output_error` — guaranteed unreachable by Phase 1.3

### 5.2 WDL transpiler

Simplify around explicit types from `OutputSpec` and normalized step bindings.

Potential cleanup:

- remove `_infer_output_type` — `OutputSpec.type` is guaranteed to be present
- remove `_dict_infer_output_type` — same reason
- less dict/object dual handling in `_binding_to_wdl_expr`
- stronger assumption that merges are already impossible

### 5.3 Nextflow transpiler

Same goals:

- remove support for non-normalized forms
- reduce dict/object dual handling
- stop using merge semantics as fallback behavior

### Outcome

Transpilers become smaller and more declarative.

---

## Phase 6 — Rationalize parser and lowering duplication

### 6.1 Workflow parser cleanup

In `syntax/wf/parser.py`:

- remove unused helpers (`_until`, `_find`)
- reduce repeated token checks by introducing small predicate helpers
- consider restructuring precedence parsing for readability

### 6.2 Lowering cleanup

In `ir/lower.py`:

- split builtin matching/desugaring from general AST lowering
- isolate chain desugaring into its own helper module
- reduce repeated block-binding lowering logic between `_lower_tree_impl()` and `lower_expr(... block ...)`
- remove dead `_record_field_names()`

### 6.3 Builtin handling

Current builtin handling is spread across parser-adjacent matchers, checker logic, and lowerer logic.

Introduce a shared builtin representation layer, e.g.:

- `syntax/wf/builtins.py` keeps AST pattern matching only
- new `semantic/wf/builtins.py` for semantic behavior
- new `ir/builtins.py` for lowering/forcing behavior

This avoids implicit coordination across modules.

---

## Phase 7 — Introduce a proper loader/cache service

Currently import/file loading and definition caching are spread across checker, lowerer, and forcer.

### Proposed service

`python/swl/loader.py` or `python/swl/project.py`

Responsibilities:

- read files
- resolve relative imports
- cache parsed tasks/workflows
- cache semantic checks
- cache lowered IR
- provide source text and provenance consistently

### Why this matters

Right now:

- checker reads files directly
- forcer reparses tasks
- lowerer caches workflow bodies separately
- multiple stages know too much about path resolution

Centralizing this would remove redundant parsing and reduce inconsistencies.

---

## Phase 8 — CLI and package surface cleanup

### `python/swl/compile.py`

Refactor into:

- a small library API module, e.g. `python/swl/api.py`
- a thin CLI wrapper

Recommended API:

- `compile_workflow(input_path, output_path=None)`
- `load_workflow(path)`
- `force_workflow(path)`
- `transpile_dag(path, target)`

If `UserError` is intended, make pipeline stages raise it consistently. Otherwise remove it and keep one error model.

### `__init__.py`

Export only stable public entry points.

---

## Recommended implementation order

1. **Delete generated artifacts and dead code**
2. **Normalize binding/DAG contract**
3. **Refactor `ir/force.py`**
4. **Refactor `semantic/wf/check.py`**
5. **Centralize type utilities**
6. **Simplify transpilers**
7. **Introduce loader/cache service**
8. **Clean CLI/public API**

This order matters because stricter contracts should come before module splitting in the transpilers.

---

## Concrete file-level plan

### Remove / delete

- `python/.pytest_cache/`
- `python/swl/**/__pycache__/`

### Re-organize

Move eval scripts to python/swl/eval/

  - `python/swl/eval.py`  -> python/swl/eval/syntax_wf.py
  - `python/swl/eval_ir.py`  -> python/swl/eval/ir.py
  - `python/swl/eval_force.py` -> python/swl/eval/force.py
  - `python/swl/eval_task.py` -> python/swl/eval/syntax_task.py
  - `python/swl/eval_task_semantic.py` -> python/swl/eval/semantic_task.py
  - `python/swl/eval_wf_semantic.py` -> python/swl/eval/semantic_wf.py

Update test.sh with new paths.

### Shrink / split

- `python/swl/semantic/wf/check.py`
- `python/swl/ir/force.py`
- `python/swl/transpile/cwl/emit.py`
- `python/swl/transpile/wdl/emit.py`
- `python/swl/transpile/nf/emit.py`
- `python/swl/ir/dag.py`
- `python/swl/ir/lower.py`

### Consolidate into single binding module

- `python/swl/ir/dag.py:_binding_to_dict`
- `python/swl/ir/dag.py:_binding_to_binding_dict`
- `python/swl/ir/dag.py:_nested_field_to_dict`
- `python/swl/ir/force.py:_binding_to_public_dict`
- `python/swl/transpile/cwl/emit.py:_canonical_binding`
- `python/swl/transpile/cwl/emit.py:_field_chain`

### Introduce

- `python/swl/ir/binding.py` — replaces all above
- `python/swl/ir/forcing/` package
- `python/swl/semantic/wf/imports.py`
- `python/swl/semantic/wf/infer.py`
- `python/swl/semantic/wf/signature.py`
- `python/swl/semantic/wf/scope.py`
- `python/swl/semantic/wf/bashvars.py`
- `python/swl/types.py` or equivalent
- `python/swl/loader.py`
- `python/swl/api.py`

---

## Testing strategy for the refactor

Before major edits, add characterization tests around current behavior.

### Priority test areas

1. **Compiler contract tests**
   - final DAG contains no `merge`
   - mapped steps always have `scatter` and `broadcast`
   - outputs always have explicit `type`
   - optionality preserved in inputs/outputs/params

2. **Binding round-trip tests**
   - binding object -> JSON -> binding object
   - same for complete DAGs

3. **Semantic inference tests**
   - record input inference
   - table input inference
   - `map` and `map_by` behavior
   - imported workflow signatures

4. **Transpiler tests**
   - transpilers consume only canonical DAGs
   - fail clearly on illegal DAGs

5. **Cleanup regression tests**
   - no reliance on removed alternate binding encodings

### Refactor safety tactic

For Phases 1–5, prefer:

- extract function/module
- add tests
- switch callers
- delete old path

rather than large rewrites in one pass.

---

## High-priority quick wins

If only a small first pass is desired, do these first:

1. delete caches and pyc files from repo
2. remove unused methods/classes
3. fix duplicate/unreachable branches
4. centralize type helpers (`_as_array_type`, optional handling, SWL type normalization)
5. stop transpilers from handling non-canonical binding variants
6. extract task/workflow definition loading from `ir/force.py`

These alone should noticeably reduce bloat.

---

## Success criteria

The refactor is successful when:

- `ir/force.py`, `semantic/wf/check.py`, and `transpile/cwl/emit.py` are substantially smaller
- final DAG JSON has one canonical binding encoding
- transpilers do not need to infer workflow output types or mapped-port semantics
- task/workflow loading is not duplicated across stages
- dead debug/cache code is removed from the package tree
- type conversion logic has a single source of truth
- new contributors can understand each stage without reading a 1000-line file

Size of codebase before:

```
$ cloc python/
github.com/AlDanial/cloc v 2.06  T=0.05 s (755.9 files/s, 147447.1 lines/s)
-------------------------------------------------------------------------------
Language                     files          blank        comment           code
-------------------------------------------------------------------------------
Python                          37           1067            140           6197
Markdown                         1              3              0              5
-------------------------------------------------------------------------------
SUM:                            38           1070            140           6202
-------------------------------------------------------------------------------
```

---

## Suggested first PR breakdown

### PR 1: hygiene + dead code removal
- remove caches/artifacts
- remove unused helpers/classes
- fix obvious duplication/unreachable code
- remove `_nested_field_to_dict` / `_field_chain` duplication

### PR 2: canonical binding module
- add normalized binding representation
- consolidate `_binding_to_dict`, `_binding_to_binding_dict`, `_binding_to_public_dict` into one
- migrate `DAG.to_dict()` / `from_dict()`
- add round-trip tests

### PR 3: forcing split
- extract finalization, merge handling, tool definition loading

### PR 4: semantic checker split
- extract imports, scope, inference, signature assembly

### PR 5: type utilities centralization
- migrate transpilers and forcer to shared type helpers

### PR 6: transpiler simplification
- remove fallback normalization logic from CWL/WDL/NF emitters (`_canonical_binding`, `_infer_output_type`, `_field_chain`, `ports` fallback, ad hoc `InputSpec`)

This sequence keeps risk manageable while steadily improving maintainability.
