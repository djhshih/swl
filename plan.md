# Remaining spec-compliance implementation plan

## Goal
Track only the concrete implementation work still remaining in `python/swl/` so the code continues to converge toward `spec.md` and `new.md`.

This plan is intentionally code-focused and excludes items already completed.

---

## Priority 1: add stronger semantic and force/DAG regression tests

### Why
The remaining work is concentrated in edge semantics where narrow regressions are easy to introduce.

### Concrete changes

#### 1. Semantic tests
Modify/add in:
- `python/swl/semantic/wf/test_check.py`

Add focused tests for:
- concrete batch root classification and column propagation from mapped callees
- remaining edge cases around imported workflows and partial application
- any newly supported or newly constrained table-update behavior
- `map_by` behavior, whether implemented or explicitly unsupported

#### 2. IR/force tests
Modify/add in:
- `python/swl/ir/test_force.py`
- `python/swl/ir/test_force_codec.py`

Add stronger tests for:
- logical mapped table-source preservation in the canonical DAG
- imported workflow mapping and partial applications inside `map`
- mapped-step schema preservation
- DAG codec round trips for mapped metadata and logical sources

These tests should stay SWL/canonical-DAG focused. Backend flattening expectations belong in CWL emitter tests.

### Required behavior after this step
- remaining semantic changes are covered by narrow regression tests before and after implementation
- canonical DAG invariants are enforced independently from CWL lowering

---

## Priority 3: review CLI/debug surfaces for the updated canonical model

### Why
The semantic and force layers now expose more explicit workflow/root typing and canonical mapped-step metadata. CLI/debug surfaces should reflect that clearly.

### Concrete changes

Review and adjust if needed:
- `python/swl/compile.py`
- `python/swl/eval_force.py`
- `python/swl/eval_ir.py`
- `python/swl/eval_wf_semantic.py`

Ensure:
- canonical DAG output is described clearly as symbolic/logical
- mapped-step source/schema information is visible enough for debugging
- validation entry points are discoverable where appropriate

### Required behavior after this step
- debug surfaces accurately reflect the current canonical representation
- users can inspect workflow typing and mapped DAG structure without confusion

---

## Priority 4: Implement table update semantics

### Why
`tab // rec`, `rec // tab`, and `tab // tab` are still not implemented. They currently fail explicitly, which is better than silent misbehavior, but this remains a major spec gap.

### Concrete changes

#### 1. Decide supported semantics
Review:
- `spec.md`
- `new.md`

Then implement or further constrain behavior for:
- `tab // rec`
- `rec // tab`
- `tab // tab`

Modify:
- `python/swl/semantic/wf/check.py`
- `python/swl/ir/force.py`
- possibly `python/swl/ir/lower.py`

If full implementation is deferred, tighten surface-level checks and errors so the unsupported boundary is completely clear.

### Required behavior after this step
- table-update behavior either works per spec or fails with precise, intentional errors
- no table-update form is silently treated like ordinary record merge when that would violate `tab` semantics

---

## Priority 5: Implement `map_by`

### Why
`map_by` is in the spec but still explicitly unsupported.

### Concrete changes

#### 1. Choose implementation vs. continued explicit rejection
Review:
- `spec.md`
- `new.md`

If implementing:
- extend semantic typing in `python/swl/semantic/wf/check.py`
- extend builtin lowering in `python/swl/ir/lower.py`
- add force/DAG and CWL behavior as needed

### Required behavior after this step
- `map_by` is implemented consistently

---


## Priority 6: continue auditing `spec.md` against implementation with focused fixes

### Why
The largest completed gaps are closed, but remaining spec mismatches should be converted into explicit tests and narrow fixes rather than broad refactors.

### Concrete changes

For each remaining mismatch found in:
- `spec.md`
- `new.md`

Do this sequence:
1. add or tighten a focused test
2. make the smallest semantic/IR/CWL change needed
3. rerun the full unittest suite

### Required behavior after this step
- remaining work proceeds incrementally and safely
- spec/design gaps become explicit tracked tests rather than implicit known issues

---

## Code-specific file checklist

### `python/swl/semantic/wf/check.py`
- implement or constrain table-update semantics
- implement or continue explicitly rejecting `map_by`
- extend focused regression coverage as behavior changes

### `python/swl/ir/lower.py`
- support `map_by` if implemented
- otherwise keep explicit unsupported handling aligned with semantic checks

### `python/swl/ir/dag.py`
- refine canonical logical table-source representation if needed
- preserve round-trip stability for mapped-step metadata

### `python/swl/ir/force.py`
- refine canonical logical table-source representation if needed
- implement or constrain table-update semantics
- preserve mapped-step logical/source/schema invariants

### `python/swl/transpile/cwl/emit.py`
- extend only if needed for any newly implemented `map_by` or table-update behavior
- keep backend flattening concerns in the emitter layer

### `python/swl/compile.py`
- verify compiled output messaging/behavior matches the canonical DAG model

### `python/swl/eval_ir.py`, `python/swl/eval_force.py`, `python/swl/eval_wf_semantic.py`
- keep debug output aligned with current workflow typing and canonical mapped representation

### Tests
- `python/swl/semantic/wf/test_check.py`
- `python/swl/ir/test_force.py`
- `python/swl/ir/test_force_codec.py`
- `python/swl/transpile/cwl/test_emit.py` when backend behavior changes

---

## Recommended implementation order

1. add stronger focused tests for remaining semantic/force gaps
2. decide and implement/constrain table-update semantics
3. decide whether to implement `map_by`
4. refine canonical table-source representation if needed
5. review CLI/debug surfaces
6. continue spec-audit work incrementally

---

## Definition of done

The remaining work is done when all of the following are true:
- canonical DAG preserves logical mapped table sources clearly and stably
- table-update semantics are either implemented per spec or explicitly constrained with precise errors
- `map_by` is either implemented consistently or rejected consistently
- regression tests cover the remaining semantic and canonical-DAG edge cases
- CLI/debug surfaces accurately reflect the current workflow typing and canonical DAG model
