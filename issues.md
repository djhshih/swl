# Issues

## Current implementation limits

- `semantic/task/type.py` now has a `Task -> TaskSignature` bridge, but the semantic model is still task-centric and not yet integrated with a richer workflow semantic model.
- Workflow semantic import resolution is now partially implemented, but only for direct `name = import "file"` bindings at the outer block.
- Current workflow semantic checking still only performs direct type compatibility checks on explicit `chain` expressions by imported task names. It does not yet perform full semantic checking for application/update-based composition.
- Workflow input inference now uses demand-driven symbolic open-record analysis instead of a hard-coded workflow input shape. This is a better direction, but it is still approximate and not yet a full semantic model.
- DAG construction and interpretation are not yet implemented.
- The workflow semantic layer now supports importing `.swl` workflows in addition to `.sh` tasks, but its workflow signatures are still approximate.
- The workflow semantic layer still has major limitations: it does not yet track precise value provenance or full record-flow semantics for arbitrary workflow expressions.
- The current checker does not yet implement the newly-settled lazy partial-application model. It still approximates some partial applications as if they were immediate applications during inference.
- The current checker does not yet implement scalar-to-record lifting for task application, where a scalar argument should be lifted into a record field named after the first declared task input.
- IR representation strategy has been discussed but not implemented yet. Current recommendation: semantic IR should be tree-shaped, while forced execution IR should be graph-shaped / DAG-like.
- We should keep a note of the tradeoff against a decorated-AST approach: a decorated AST is simpler initially, but imports, chain normalization, closures, lazy application, and later forcing introduce semantic objects that differ enough from raw syntax to justify a small lowering/normalization step.
- We should continue to validate whether a single semantic `IRImport(kind=...)` node is sufficient, or whether task/workflow imports eventually need separate IR node classes. Current recommendation is to use a single import node in semantic IR and branch on kind only during forcing.
- The current checker does not yet enforce that a workflow must evaluate to a function.
- Workflow output inference does not yet implement the newly-specified chain-output rule: when a workflow evaluates to a chain, outputs should be the union of the output variables from left to right.
- Task semantics do not yet enforce that all task output params must have defaults.
- Task output defaults may include glob patterns, but glob-aware validation/representation is not yet implemented explicitly.
- `bash.py` is intended to support pre-runtime and runtime bash validation after interpolation, but this two-stage validation model is not yet implemented.
- All issues are intended to be errors unless explicitly specified otherwise, but the current code still distinguishes some categories operationally (`chain_errors` vs `issues`).
