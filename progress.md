# Progress

- Implement table update semantics (`tab // rec`, `rec // tab`, `tab // tab`). Checker has `_merge_table_record` but force.py rejects table updates at DAG construction time; test is commented out.
- Allow binding statements in a block to reference each other, but no recursion. Spec was updated but there is no forward-reference validation — the evaluator would fail silently at name resolution.
- Audit remaining `new.md` / `spec.md` gaps and convert each into focused tests plus narrow fixes.
- Add stronger test coverage for concrete batch root typing, table-column propagation, and mapped table source round-trip stability.

### Remaining gaps (not fully guaranteed)

| Gap | Status |
|-----|--------|
| Final-DAG merge freedom | `_flatten_merge_value` in force.py handles most cases; merge bindings should not survive to DAG output but edge cases may remain. |
| Full optionality propagation | `Input` has `optional` field and DAG serializes it, but some paths may miss the `?` from task annotation parsing. |
| Universal mapped-port classification | `_mapped_step_bindings` sets `scatter`/`broadcast` for every mapped step; verify every input lands in exactly one of the two. |
| Non-saturating record bindings | Records that do not directly feed a task call interface may still reach the DAG as `record` bindings. CWL rejects these. |
| Nested field projection (`Field(Field(...))`) | CWL transpiler still rejects these; WDL and NF handle them inline. |
| Record merge in workflow outputs | Merge flattening should cover this but untested for all merge-tree shapes at the output level. |
| Output descriptions on workflow outputs | `OutputSpec.desc` is now supported in the DAG model, but current forcing populates it as `None` because workflow-output description inference is not yet implemented. |
| Full optionality propagation beyond basic `?` typing | Better than before: `OutputSpec.optional` now exists and transpilers consume serialized output optionality. Some deeper inference paths may still need auditing. |
| Output type inference for complex constructed outputs | Basic output typing is now compiler-owned. Deeply nested/constructed record-field outputs still fall back to `None` in force-time inference and may need future refinement if such outputs become part of the final DAG contract. |
