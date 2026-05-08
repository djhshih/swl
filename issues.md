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
- IR representation strategy is now partially implemented: semantic IR is tree-shaped, but forced execution IR / DAG construction is still missing.
- We should keep a note of the tradeoff against a decorated-AST approach: a decorated AST is simpler initially, but imports, chain normalization, closures, lazy application, and later forcing introduce semantic objects that differ enough from raw syntax to justify a small lowering/normalization step.
- The semantic IR no longer uses an import-specific node; imports are reduced eagerly to cached `Function` values. Remaining question: whether `Closure` should also be folded into a more unified function-like representation.
- Workflow output inference does not yet implement the newly-specified chain-output rule: when a workflow evaluates to a chain, outputs should be the union of the output variables from left to right.
- Task output defaults may include glob patterns, but glob-aware validation/representation is not yet implemented explicitly.
- `bash.py` is intended to support pre-runtime and runtime bash validation after interpolation, but this two-stage validation model is not yet implemented.
- All issues are intended to be errors unless explicitly specified otherwise, but the current code still distinguishes some categories operationally (`chain_errors` vs `issues` / `errors`).
- Scope validation is now done semantically. The language rule is that nested scopes may shadow outer scopes, while duplicate bindings in the same immediate scope are errors.
- The current checker's lazy partial-application model is still too field-set/signature oriented: closures remember mostly bound field names rather than bound values/provenance, and saturated applications become approximate computation summaries rather than structured lazy application values.
- The compiled DAG JSON is now executor-oriented and self-contained enough to avoid rereading `.swl` / `.sh` files at execution time, but some workflow-level input metadata is still lossy. In particular, when a root workflow reaches forcing through an inferred lambda/chain shape rather than a directly imported task signature, external input entries may lack precise type/description metadata.
- The compiled task payload now intentionally omits raw annotation section structure in favor of normalized `inputs` / `outputs` / `run`, but there is still no frozen executor contract for how interpolation expressions in defaults and run parameters should be evaluated against runtime bindings.
- `pipe.swl` and `function.swl` are intended to be equivalent, but compilation currently diverges because `pipe.swl` reaches forcing as an explicit `ir.Chain` while `function.swl` reaches forcing as nested `Apply` + `Update` structure ending in a `Merge`. Today `force.py` canonicalizes chain outputs better than explicit record-update outputs, so `function.swl` still serializes its final outputs under a nested `result` merge instead of the same flat output map used by `pipe.swl`. The first fix should be in `ir/force.py`: flatten final record merges into output maps and preferably canonicalize record merges eagerly. A later optional cleanup is to teach `ir/lower.py` to recognize explicit workflow-composition patterns and lower them toward the same canonical chain-oriented shape.
