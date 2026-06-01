# Remaining spec-compliance implementation plan

## Goal
Track only the concrete implementation work still remaining in `python/swl/` so the code continues to converge toward `spec.md` and `new.md`.

This plan is intentionally code-focused and excludes items already completed.

---

## Priority 5: Implement CWL transpilation for `map_by`

### Why
`map_by` is now implemented in semantic checking, lowering, and forcing, and canonical DAGs preserve grouped mapping via mapped-step metadata (`map.group_by`). The remaining gap is backend lowering: `python/swl/transpile/cwl/emit.py` still rejects grouped mapped steps explicitly.

### Scope for this plan
This section is only about CWL transpilation for workflows involving `map_by`.
It does **not** reopen `map_by` language semantics unless emitter work exposes a concrete spec gap.

### Constraints from the current design
- Canonical DAG must remain logical and preserve grouped mapping metadata.
- Grouping semantics stay in the canonical model as:
  - mapped execution
  - logical table source
  - `map.group_by`
- Backend-specific flattening must stay in CWL emission.
- If a faithful CWL lowering cannot be produced for a case, emission must fail with an explicit intentional error rather than silently degrading to ordinary `map`.

### Concrete changes

#### 1. Decide the backend strategy for grouped mapping in CWL
Review:
- `spec.md` (`map_by` section)
- `new.md`
- current `python/swl/transpile/cwl/emit.py`
- current canonical mapped-step shape in `python/swl/ir/dag.py`

Pick one explicit strategy and encode it in tests before broad implementation:

Option A: grouped-subworkflow lowering in the emitter
- emit a preprocessing/grouping stage that constructs per-group inputs
- emit the mapped grouped callee over those grouped inputs
- preserve one output row per unique key
- preserve the grouping key in emitted outputs

Option B: continue rejecting some shapes explicitly, but support a narrow faithful subset
- for example, only support root/grouped mapped workflow steps whose grouped input/output schemas are statically known enough to synthesize CWL wiring safely
- reject all other grouped shapes with precise errors

This choice should be recorded in comments in `python/swl/transpile/cwl/emit.py` and reflected in tests.

#### 2. Make grouped mapped-step metadata first-class in emitter validation
Modify:
- `python/swl/transpile/cwl/emit.py`

Current behavior:
- grouped mapped steps are rejected immediately

Change validation so it distinguishes:
- supported grouped-map shapes
- intentionally unsupported grouped-map shapes

Required checks should include:
- `map.group_by` exists and is a string
- grouping key exists in the mapped step input schema / logical source schema
- grouping key is preserved in the mapped step output schema
- required logical table-source metadata is present for emitter-side lowering

If any of those are missing, raise a precise `ValueError` that explains which invariant was violated.

#### 3. Introduce an emitter-side representation for grouped table preparation
Modify as needed:
- `python/swl/transpile/cwl/emit.py`
- possibly helper structures local to CWL emission only

The emitter likely needs an internal helper that can derive, from a mapped step:
- the grouping key name
- the logical input table columns
- the grouped per-key input object shape
- the expected per-group output object shape

Keep this representation emitter-local.
Do not move grouping flattening back into forcing or canonical DAG normalization.

#### 4. Define how grouped inputs are represented in emitted CWL
Implement a concrete emitted shape for grouped inputs.
At minimum, document and encode how CWL will represent:
- the set of unique grouping-key values
- each grouped slice passed to the mapped callee
- non-key grouped columns as arrays aligned within the group
- the grouping key itself in each grouped invocation

The emitted shape must reflect the spec requirement that grouped slices are `rec`, not `tab`, even if some grouped fields are arrays.

That means the CWL lowering should conceptually feed each grouped invocation an object like:
- `{ key: scalar, col1: [..], col2: [..], ... }`
rather than a backend-level table abstraction.

#### 5. Lower grouped mapping to CWL workflow structure
Modify:
- `python/swl/transpile/cwl/emit.py`

Implement emitter logic that turns one canonical mapped step with `group_by` into an emitted CWL workflow structure.
The exact structure depends on the chosen strategy, but should include the equivalent of:
1. derive grouped inputs from the incoming logical table columns
2. invoke the grouped callee once per group
3. collect outputs back into arrays / row-like columns representing one output row per unique key

The grouped lowering must preserve:
- one emitted output element per unique key value
- grouping-key preservation
- stable scatter/wiring for grouped outputs
- compatibility with downstream field access like `calls.bcf`

#### 6. Decide how to implement the grouping operation itself in CWL
Review whether current emitter patterns can express grouping directly or whether a helper tool/subworkflow is needed.

Possible implementation directions:
- synthesize a CWL `ExpressionTool` that groups parallel input arrays by key and emits grouped records
- synthesize a helper `Workflow`/tool layer specifically for grouping
- encode grouping using existing CWL expression-capable nodes if already present in the project’s emission model

Whichever path is chosen, keep the implementation localized to emission.
Do not change canonical DAG semantics just to fit CWL.

#### 7. Preserve output cardinality and key retention explicitly
Add explicit emitter checks and output shaping so that:
- emitted grouped results have one element per unique key
- the grouping key is included in emitted grouped outputs
- downstream consumers can access grouped result columns with ordinary field access semantics from the canonical model

If the grouped callee output schema omits the key, emission should fail explicitly even if semantic checking already catches it. The emitter should defend its own invariants.

#### 8. Keep partial-root `map_by` workflows working through CWL emission
Current code now supports root partial `map_by` in canonical DAGs.
Ensure emitter behavior covers canonical DAGs of the form:
- root workflow is effectively `map_by f "k"`
- mapped step source is a root input like `{'source': 'input', 'name': 'xs'}`
- `map.group_by` survives codec round trips and is consumed by emission

This is a distinct regression surface from in-body grouped mapping and should have its own tests.

### Required test plan

#### `python/swl/transpile/cwl/test_emit.py`
Add focused tests for:
- grouped mapped step currently supported shape emits successfully
- root partial `map_by` emits successfully once supported
- grouped mapped workflow emits grouping helper/preprocessing structure as expected
- grouping key must exist in input schema for emission
- grouping key must be preserved in output schema for emission
- unsupported grouped shapes fail with precise errors
- grouped outputs remain consumable by downstream field access in emitted CWL

#### `python/swl/ir/test_force.py`
Add or tighten tests ensuring forced DAGs expose the metadata the emitter relies on:
- `map.group_by`
- logical mapped source
- `input_schema`
- `output_schema`

#### `python/swl/ir/test_force_codec.py`
Keep round-trip coverage for:
- partial root `map_by`
- in-body grouped mapped steps
- preserved `group_by` metadata after JSON round trip

### Required behavior after this step
- CWL emission supports the chosen faithful subset of `map_by`
- any unsupported grouped-map shape fails with a precise, intentional error
- canonical DAG semantics remain logical and unchanged
- grouped mapping is not silently degraded into ordinary scatter over rows

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

