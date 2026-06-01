# Issues

Issues are concrete gaps between the spec and the implementation, or bugs that violate the "fail fast" principle. Design choices and vague observations have been removed.

---

## Checker and semantic analysis gaps

### Import resolution only works at the top-level block

**Spec:** `import: str -> fun` should be usable anywhere an expression is valid.

**Code:** `semantic/wf/check.py:_load_imports()` (line 172) only scans the outer block for `name = import("path")` patterns. An `import` call inside a lambda body, nested block, or as a subexpression is parsed as a regular function application but never resolved by the semantic checker or the lowerer. Such imports will silently fail at forcing time when the unresolved `Function` node is encountered.

### Type compatibility is only checked for explicit `|` chain syntax

**Spec Sec 2 (Type Compatibility):** Both pipeline (`a | b`) and record union (`r1 // r2`) must be type-checked. For pipelines: matching output/input field types must be compatible per the type matrix. For record unions: same-named fields must have compatible types.

**Code:** `Checker._check_chains()` (line 246) only walks the AST for `NodeType.chain` nodes and calls `TypeChecker.check_chain()`. No type checking is performed for:
- `Apply`-based composition (`f x` where `x` is the result of another call)
- `Update` expressions (`r1 // r2`) — `_merge_update_values` merges without checking type compatibility
- Chain desugaring in the lowerer (`_chain_to_lambda_block`) doesn't re-check types since it operates after the checker

### Partial application model drops provenance information

**Spec (Function Application):** "Partial application returns a new function." The partial function should preserve what's been bound.

**Code:** `_application_result()` in `check.py` (line 860) creates `ClosureValue` instances that remember bound field names but not the actual values or their provenance. For example, if `task1` has inputs `{a, b, c}` and you partially apply with `{a: some_record.field}`, the resulting closure records `{'a': UnknownValue}` rather than `{'a': Field(some_record, 'field')}`. This means later compilation stages (lowering → forcing) may produce different step call bindings than the semantic checker assumed. The checker also eagerly creates `ComputationValue` (treating a function as applied) even when the bound record is an `OpenRecord` rather than a concrete `ClosedRecord`, which makes inference approximate.

### DAG no circularity check

**Spec Sec 3:** "DAG graph will be constructed and checked for circularity."

**Code:** `ir/dag.py:DAG` stores `steps` with dependency lists but has no cycle detection. `ir/force.py` can produce a DAG with circular step dependencies if a workflow expresses one (e.g., through mutually-recursive imports or self-referential bindings). The only circularity check is for file imports (`Checker._loading` stack), not for step-level dataflow.

---

## Execution and serialization gaps

### DAG JSON has no contract for evaluating interpolation defaults

**Spec (Type Annotations, Features):** String interpolation in defaults (`${outbase}.bam`) must be evaluated at execution time.

**Code:** `force.py:_build_task_param()` stores parsed default values as interpolation AST nodes (`Word` → `Literal`/`Var`/`Expr`) serialized to dicts (line 842-851: `_interp_to_dict`). The executor receives these AST fragments but there is no documented contract for how to evaluate them against runtime bindings. The executor must know to substitute `Var` references with input/output values and evaluate `Expr` fragments, but no executor interface exists.

### Workflow-level input metadata is lossy in compiled DAG

**Spec Sec 4 (Inference):** Inferred inputs should carry proper type and description metadata.

**Code:** When the root workflow reaches forcing through an inferred shape (lambda/chain) rather than a directly imported task signature, `force.py:_refine_input_metadata()` (line 556) attempts to merge type info from step input specs, but if the top-level input doesn't match any task input name, the `Input` is created by `_input()` (line 723) with `type=None` and `desc=None`. This means inferred workflow inputs that are only used as lambda parameters or record fields lack type/desc in the final DAG.

### `pipe.swl` and `function.swl` produce different DAGs despite semantic equivalence

**Spec (Pipeline):** `A | B | C` desugars to `\x -> a = A x; b = B (x // a); c = C (x // a // b); a // b // c`.

**Code:** `lower.py:_chain_to_lambda_block` correctly desugars chains. However, the equivalent hand-written form using nested `Apply` + `Update` reaches `force.py` as a different IR shape. `force.py:_final_outputs()` handles the chain case's flat merge better than the equivalent update-tree case. The result is that `function.swl` serializes outputs under a `result` merge node instead of the same flat output map used by `pipe.swl`. The root cause is in `_normalize_output_value()` (line 793) which doesn't fully flatten recursive `Merge` structures.

### No glob validation for task output defaults

**Spec (Type Annotations):** Output parameter defaults can include file paths.

**Code:** Task output defaults like `bam/*.bam` are parsed as interpolation words but never validated for glob correctness or syntax. If a glob pattern is malformed, the error only surfaces at task execution time (if at all) rather than at compile time.

### `bash.py` two-stage validation not integrated

**Spec (Compile-time Checks):** Catch as many issues as possible statically.

**Code:** `syntax/task/bash.py` parses bash scripts into `Assignment`/`Command` structures with interpolation analysis, but no part of the pipeline uses it for validation. The compiled DAG stores the raw bash body as `script`, and no static analysis checks for undefined variable references, missing tool invocations, or shell syntax errors before execution.

---

## Code structure and correctness issues

### Parser doesn't allow lambda on the right side of `|`

**Spec:** `chain ::= expr ( "|" expr)+`. `expr ::= ... | lambda`. Precedence: `->` > `|`. So `a | \x -> body` should parse as `a | (\x -> body)`.

**Code:** `_parse_chain_expr()` (line 185) calls `_parse_update_expr()` for the right operand, which eventually reaches `_parse_term()`. `_parse_term` only handles `str`, `num`, `id`, `lbrace`, `lparen` — not `bslash` (`\`). So any lambda on the right side of a `|` operator causes a parse error. The lambda must be parenthesized (`a | (\x -> body)`) to parse correctly.

### Parser doesn't enforce spec's block grammar

**Spec:** `block ::= (binding eol)* expr` — only bindings (with trailing eol) before the final expression. Non-binding expressions in the middle of a block are invalid.

**Code:** `_parse_block()` collects all expressions separated by `eol` tokens, only validating that the final expression is not a `Binding`. The input `a = 1\nb\nc = 2\nd` parses without error even though `b` (a bare identifier expression) appears in the middle of the block, which violates the spec grammar.

### Parser is more permissive than spec for `apply` and `get` left-hand sides

**Spec:** `apply ::= name ws expr` — the function position must be a bare `name`. `get ::= name ("." name)+` — the record position must be a bare `name`.

**Code:** `_parse_apply_expr()` (line 221) allows `fun`, `id`, `get`, `apply`, `chain` on the left side. `_parse_get_expr()` (line 245) also allows `rec` on the left side. While this is likely an intentional ergonomic choice, it means the implementation accepts programs the spec grammar rejects.

### `syntax/wf/node.py`: `Expr.__repr__` accesses subclass-only attributes on the base class

`__repr__` switches on `self.type` then accesses `self.name`, `self.value`, `self.body`, `self.id`, `self.rec`, `self.member`, `self.left`, `self.right`, `self.fun`, `self.arg`, `self.param` — all defined only on specific subclasses. Adding a new `NodeType` without updating `__repr__` causes `AttributeError` at display time. A dispatcher method per node class would be more maintainable.

### `ir/node.py`: `Unknown` sentinel can propagate to forcing

`ir/lower.py` returns `ir.Unknown()` (line 151) as the fallback for unhandled AST node types. No handler exists in `force.py:force_value()` for `ir.Unknown`, so it raises a generic `ValueError` at forcing time. The error doesn't point back to the original AST location or node type that produced it.

### `ir/force.py`: Unresolved variable names silently become workflow inputs

`force_value()` for `ir.Name` nodes (line 93) checks the `ForceEnv`; if the name is absent, it creates a new `Input` without warning. A typo in a variable name (e.g., `algn_results` instead of `align_results`) silently becomes a workflow input rather than a compile error. This violates the "fail fast" principle.

### `ir/force.py`: Run parameters with missing values silently dropped

`_normalize_task_run()` (line 459) only includes run parameters that can be projected from the bound record. A run parameter like `cpu` with no default that isn't in the bound record is silently omitted. The DAG compiles successfully, but the executor will receive no value for the required run parameter.

### `ir/dag.py`: `StepCall` and `MappedStep` serialization type confusion

`to_dict()` uses `getattr(step, 'map', None)` on all steps, but `map`, `input_schema`, `output_schema` only exist on `MappedStep`. `from_dict` selects the class by checking for the `map` key — if a `StepCall` somehow acquires map data, it's silently dropped. The dual-class design is fragile because serialization code must know about both types.

### `ir/dag.py`: `_run_value_from_dict` called twice per item

The dict comprehension in `from_dict` (line 152-155) calls `_run_value_from_dict` in both the value expression and the `if` filter, doubling the work. The function also returns `None` for multiple semantically different reasons (param not found, value equals default, unparseable value), making the contract unclear.

### Duplicate pattern-matching for `import`, `map`, and `map_by`

Both `semantic/wf/check.py` and `ir/lower.py` independently implement pattern-matching for `import("path")`, `map f xs`, and `map_by f key xs` using nearly identical AST-walking logic. Any change to the AST representation or built-in detection rules must be mirrored in both files.

### `semantic/task/type.py`: `Param.__init__` type annotations mismatch

`Param.__init__` declares `typ: TypeKind`, `default: str`, and `desc: str` as non-optional, but `None` is regularly passed for all three throughout the codebase (check.py lines 423, 431, 523; lower.py lines 384, 386, 401). This produces LSP type errors in every file that constructs `Param` with missing metadata.

### `compile.py`: All exceptions handled identically

`compile_workflow` propagates exceptions, but the `__main__` block catches everything and prints to stdout with `os.EX_DATAERR`. There's no `--verbose` flag, no distinction between user errors (bad input) and internal bugs, and tracebacks go to stdout instead of stderr.

### `WorkflowCheck.chain_errors` and `WorkflowCheck.issues` are vestigial

Both properties in `check.py:WorkflowCheck` (line 118-123) return `self.errors`. They are unused aliases from an earlier design that distinguished error categories. They should be removed to avoid confusion about whether different validation paths exist.
