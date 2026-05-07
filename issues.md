# Issues

## Open semantic questions

- Interpolation values are currently preserved as syntax nodes (`Word`, `Var`, `Expr`) in task defaults. This seems correct, but the exact resolution phase still needs to be defined clearly.
- It is not yet fully settled whether task output params must always have defaults or whether they may be derived only from body/runtime behavior.
- `bash.py` exists as an optional analyzer, but current semantics should probably treat annotation metadata as authoritative and body text as opaque.

## Current implementation limits

- `semantic/task/type.py` now has a `Task -> TaskSignature` bridge, but the semantic model is still task-centric and not yet integrated with a richer workflow semantic model.
- Workflow semantic import resolution is now partially implemented, but only for direct `name = import "file"` bindings at the outer block.
- Current workflow semantic checking still only performs direct type compatibility checks on explicit `chain` expressions by imported task names. It does not yet perform full semantic checking for application/update-based composition.
- Workflow input inference now uses demand-driven symbolic open-record analysis instead of a hard-coded workflow input shape. This is a better direction, but it is still approximate and not yet a full semantic model.
- DAG construction and interpretation are not yet implemented.
- The workflow semantic layer now supports importing `.swl` workflows in addition to `.sh` tasks, but its workflow signatures are still approximate.
- The workflow semantic layer still has major limitations: it does not yet track precise value provenance, robust partial application semantics, or full record-flow semantics for arbitrary workflow expressions.
