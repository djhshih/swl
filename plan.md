# Implementation Plan

Completed implementation notes have been moved to `done.md`.
This file tracks only the remaining work and the next implementation steps.

## Current state

Already implemented in `python/`:
- workflow parser
- task syntax parser/interpolation/bash analysis split
- task semantic typing and duplicate checks
- explicit semantic validation/normalization for `run` params (`cpu`, `memory`, `time`, `image`)
- workflow semantic checking with imports, scope checks, and basic lazy/partial-application approximation
- semantic IR lowering under `python/swl/ir/`
- canonical lowering to `Lambda(Block(...))`
- force-to-DAG compilation
- DAG JSON codec
- compiled task `run` metadata using normalized `type` + `value`
- diagnostics and test integration

The remaining work is now mostly about semantic precision, cleanup, and executor-oriented polish.

## Phase 1: tighten workflow semantics

### Goal
Make workflow semantic checking less approximate around lazy application, closures, and record provenance.

### Remaining issues
- closure semantics still track field availability more than concrete bound values/provenance
- computation values are still signature-oriented rather than execution-oriented
- some semantic conclusions are still inferred from field sets instead of preserved value flow
- error reporting still exposes some legacy naming (`chain_errors`) instead of one clean error surface

### Planned files
- modify `python/swl/semantic/wf/check.py`
- modify `python/swl/semantic/wf/test_check.py`
- optionally update `python/swl/eval_wf_semantic.py`

### Next steps
1. Make closure modeling more provenance-aware.
   - preserve which values/records were actually bound, not just which field names are available

2. Tighten partial-application semantics.
   - keep distinguishing function/closure values from computation/record-like values
   - continue rejecting field access on function values early
   - reduce remaining approximation gaps for nested applications

3. Simplify semantic error reporting.
   - move toward one `errors` surface
   - keep compatibility shims only if still needed by tests/tools

## Phase 2: clean up semantic IR

### Goal
Keep semantic IR minimal, canonical, and closely aligned with what forcing expects.

### Remaining issues
- `ir.Chain` still exists as a transient node even though forcing rejects it
- some older plan assumptions in code/comments still refer to superseded IR shapes
- `force.py` is still structurally large and mixes multiple concerns

### Planned files
- modify `python/swl/ir/node.py`
- modify `python/swl/ir/lower.py`
- modify `python/swl/ir/force.py`
- modify `python/swl/ir/test_lower.py`
- modify `python/swl/ir/test_force.py`

### Next steps
1. Remove transient `ir.Chain` completely.
   - if lowering always normalizes it away, the node may no longer need to exist

2. Refactor `python/swl/ir/force.py`.
   - separate forcing logic, task-definition shaping, normalization helpers, dependency extraction, and codec-related helpers more clearly
   - keep behavior unchanged while improving readability and maintainability

3. Tighten invariants around canonical IR.
   - keep failing loudly on unsupported/non-normalized forms
   - add direct regression tests where useful

## Phase 3: improve executor-facing compiled metadata

### Goal
Make compiled DAG JSON more directly usable by an executor, with less need to reinterpret raw syntax.

### Remaining issues
- task `run` metadata is normalized, but `inputs` / `outputs` may still expose more syntax-oriented structure than necessary
- compiled metadata shape should be reviewed for consistency across inputs, outputs, and run params
- current executor-facing JSON behavior should be covered more explicitly by tests

### Planned files
- modify `python/swl/ir/force.py`
- modify `python/swl/ir/dag.py`
- modify `python/swl/ir/test_force.py`
- modify `python/swl/ir/test_force_codec.py`

### Next steps
1. Review whether task `inputs` and `outputs` should carry more normalized semantic metadata.
   - keep enough information for execution
   - avoid unnecessary raw-syntax leakage where possible

2. Add regression tests for compiled metadata shape.
   - especially `run[*].type`
   - and `run[*].value`

3. Keep the compiled artifact self-contained.
   - executors should not need to reread original `.swl` / `.sh` files

## Phase 4: optional task/bash follow-up

### Goal
Decide whether to deepen task body validation beyond the current conservative syntax pass.

### Remaining issues
- bash-body analysis is intentionally conservative
- two-stage validation (pre-runtime vs runtime after interpolation) is still mostly a design direction, not a completed subsystem

### Planned files
- modify `python/swl/syntax/task/bash.py`
- optionally add new semantic/runtime validation helpers under `python/swl/semantic/task/`
- modify `python/swl/syntax/task/test_bash.py`

### Next steps
1. Decide whether pre-runtime bash validation should become a real semantic/compiler phase.
2. If yes, validate concrete task-call environments after interpolation becomes sufficiently known.
3. Keep this work separate from the core workflow/IR cleanup unless it becomes blocking.

## Recommended immediate next step

Work on **Phase 1: tighten workflow semantics**.

Why this is next:
- parser, lowering, forcing, codec, and normalized run metadata are already in place
- the biggest remaining source of incorrectness is semantic approximation around closures and partial application
- improving that layer should make both diagnostics and later executor behavior more trustworthy
