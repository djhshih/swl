# Issues

Issues are ordered by priority. **Priority 1** items are architectural — they stem from the three-pass pipeline design (checker → lowerer → forcer) where analysis results computed by one stage are discarded before the next. **Priority 2** items are missing compile-time checks that violate the "fail fast" principle. **Priority 3** items are execution/serialization contract gaps. **Priority 4** items are parser correctness bugs. **Priority 5** items are code quality.

---

## Priority 1: Architectural — three-pass pipeline has no constraint bridge

The pipeline evaluates the program three times with different value domains:

| Pass | Input | Domain | Output |
|------|-------|--------|--------|
| Checker | WF AST | `OpenRecord`, `FunctionValue`, `ClosureValue`, `ComputationValue`, ... | `WorkflowCheck` (signature + type) |
| Lowerer | WF AST | `ir.Name`, `ir.Ref`, `ir.Variable`, `ir.Apply`, `ir.Lambda` | IR tree |
| Forcer | IR tree | `dag.Input`, `dag.Record`, `dag.Field`, `dag.StepCall`, `ForcedFunction` | `DAG` |

The Checker's analysis is not embedded in the IR. The Forcer re-derives everything from scratch. Three architectural issues follow from this.

---

### P1.1: Partial application model drops provenance between checker and forcer

**Root cause:** Checker and forcer each evaluate partial application independently with different value domains. The checker produces `ClosureValue(function, bound_fields)` where `bound_fields` is `ClosedRecord` of `UnknownValue` entries — it remembers that fields are bound but not what they are bound to. The forcer produces `ForcedFunction(node, bound_value, signature)` with concrete bound values. There is no validation that the forcer's concrete evaluation agrees with the checker's abstract prediction.

**Spec violation:** "Partial application returns a new function." The checker's `ClosureValue` loses the bound values, so later stages may produce different step call bindings than the checker assumed. The checker also eagerly creates `ComputationValue` (treating a function as applied) when the bound record is an `OpenRecord`, making inference approximate.

**Proposed solution:** Thread the checker's type/constraint map through the IR so the forcer can use it:
1. Add a type-constraint field to `ir.Lambda` and `ir.Function` nodes that records which inputs the checker determined as demanded vs satisfied.
2. In `force.py:_input()`, consult this constraint map to set `type` and `desc` on newly created `Input` nodes instead of defaulting to `None`.
3. Add `force.py:_validate_bindings()` that cross-checks each step call's bindings against the checker's predicted input set and reports mismatches as compile errors.

---

### P1.2: Checker-inferred input metadata is lost between checker and DAG output

**Root cause:** The checker infers `inputs = {a: file, b: str}` and produces a `TaskSignature`. The lowerer receives this signature and attaches it to the root `ir.Lambda`, but the forcer's `_input()` method (`force.py:723`) only receives a `spec` parameter when the caller happens to pass one. The `_refine_input_metadata()` recovery pass (`force.py:556`) scans step-level task input specs for matching names — input names that are only used as lambda parameters or record fields never get type/desc metadata.

**Spec violation (Sec 4):** Inferred inputs should carry proper type and description metadata. The spec requires that "field requirements are inferred from unsatisfied task inputs and field accesses" — the checker does this inference but the result is not propagated to the DAG.

**Proposed solution:** Two options, one minimal and one thorough:

*Minimal:* Attach the checker's full `TaskSignature` dict to the root IR node so the forcer can look up each input name's type/desc spec. In `force.py:_input()`, if the name exists in the root signature, use its type/desc instead of `None`. This is a data-flow fix: the information exists in the `WorkflowCheck` but is not passed into the forcer.

*Thorough:* Add `ir.Input(name, type, desc)` as a first-class IR node. The lowerer creates `ir.Input` for names that the checker identified as demanded workflow inputs (vs names that are lambda parameters or local variables). The forcer then never creates `Input` implicitly — it only receives them from the IR. A name that was never declared as an `ir.Input` and is not in the `ForceEnv` becomes a compile error, not a silent input.

---

### P1.3: IR has no way to distinguish "workflow input" from "undefined variable"

**Root cause:** The IR has no `ir.Input` node. The lowerer produces `ir.Name(expr.name)` for any identifier not found in its environment — whether that identifier is a workflow input parameter or a typo. The forcer can't tell the difference. When it sees `ir.Name("align_inpt")` not in the `ForceEnv`, its only option is to create a new `Input`, silently accepting a typo as a workflow input.

**Spec violation:** Violates the "fail fast" principle. A misspelled variable name should produce a compile-time "undefined variable" error.

**Proposed solution (same thorough option as P1.2):**
1. Add `ir.Input(name, type, desc)` as an IR node.
2. The lowerer creates `ir.Input` for names that the checker identified as demanded workflow inputs. For lambda parameters it creates `ir.Name` (which gets bound in the `ForceEnv`). For truly undefined names (not in env, not in imports, not demanded by checker), it raises a compile error.
3. The forcer's `force_value()` for `ir.Input` returns the corresponding `dag.Input` from its inputs dict. For `ir.Name` it only looks in the `ForceEnv` — if absent, it's a forcing-time error.
4. This eliminates the fallback path where `ir.Name` becomes an implicit input.

---

### P1.4: Duplicate pattern-matching for imports and builtins across checker and lowerer

**Root cause:** Both `semantic/wf/check.py` and `ir/lower.py` independently implement pattern-matching for `import("path")`, `map f xs`, and `map_by f key xs`. The checker needs it for type inference; the lowerer needs it for IR construction. Any change to the AST structure or built-in detection rules must be mirrored in both files.

**Proposed solution:** Extract a shared pattern-matching module, e.g. `syntax/wf/builtins.py`, with functions like `match_import(expr)`, `match_map(expr)`, `match_map_by(expr)`. Both checker and lowerer import from this single source. If the detection logic later depends on checker state (e.g., resolved import tables), make these methods on a shared base class or pass the checker state as a parameter.

---

## Priority 2: Missing compile-time checks

### P2.1: No DAG circularity check

**Spec Sec 3:** "DAG graph will be constructed and checked for circularity."

**Code:** `DAG` stores `steps` with dependency lists (`deps`) but has no cycle detection. `force.py` can produce a DAG with circular step dependencies (e.g., through self-referential bindings or mutually-recursive imports not caught by the import circularity check).

**Proposed solution:** Add `DAG.validate()` that runs `topological_sort()` on the step dependency graph and raises `ValueError` if a cycle is found. Call it in `force.py:_finalize_dag()` before returning the DAG.

---

### P2.2: Run parameters with missing values silently dropped

**Spec:** Run parameters (`cpu`, `memory`, `time`, `image`) are required for execution. If they have no default, the workflow input must provide them.

**Code:** `_normalize_task_run()` (`force.py:459`) projects run param names from the bound record and silently omits any name where `_project_field` returns `_SENTINEL` (i.e., the value is not available). The DAG compiles successfully but the executor receives no value for the required run parameter.

**Proposed solution:** After the projection loop, check that every run parameter without a default was successfully projected. If any are missing, raise `ValueError` with the parameter name and the task path. Allow omission only for run parameters that have a `parsed_default` in their signature.

---

### P2.3: Type compatibility only checked for explicit `|` chain syntax

**Spec Sec 2:** Both pipeline (`a | b`) and record union (`r1 // r2`) must be type-checked.

**Code:** `Checker._check_chains()` only walks `NodeType.chain` nodes. No type checking is performed for `Apply`-based composition or `Update` expressions. `_merge_update_values` merges records without checking type compatibility of same-named fields.

**Proposed solution:**
1. In `Checker._eval_expr` for `NodeType.apply`, call `TypeChecker.check_chain` (or equivalent) on the function and argument when both resolve to known task signatures.
2. In `_merge_update_values`, when both sides are `ClosedRecord` or `OpenRecord`, validate that overlapping field names have compatible types per the spec's type matrix.
3. Add chain-desugaring validation in the lowerer: after `_chain_to_lambda_block`, walk the composed bindings and check each step's output/next-step-input type compatibility.

---

### P2.4: Import resolution only works at the top-level block

**Spec:** `import: str -> fun` should be usable anywhere an expression is valid.

**Code:** `_load_imports()` only scans the outer block for `name = import("path")` patterns. An `import` inside a lambda or nested block is parsed but never resolved. It fails at forcing time with an unresolved `Function` node.

**Proposed solution:** In `Checker._eval_expr` for `NodeType.apply`, when the function is `Identifier("import")` and the argument is a string literal, resolve the import inline. Cache the resolved `FunctionValue` as with top-level imports. In `Lowerer.lower_expr` for the same pattern, produce an `ir.Function` node directly. This makes import resolution work in any expression position.

---

### P2.5: No glob validation for task output defaults

**Spec:** Output parameter defaults are string expressions evaluated at execution time. Glob patterns are common (e.g., `bam/*.bam`).

**Code:** Output defaults are parsed as interpolation words but never validated. A malformed glob pattern (e.g., `[invalid`) is silently accepted.

**Proposed solution:** In `force.py:_build_task_param()`, when the parameter is an output with a default, apply a glob syntax check. Use Python's `glob` module or a simple regex to validate that the pattern is well-formed. Report a compile-time error for malformed patterns. This can optionally be extended in `bash.py` to verify that glob patterns reference existing files at runtime.

---

### P2.6: `bash.py` two-stage validation not integrated

**Spec (Compile-time Checks):** Catch as many issues as possible statically.

**Code:** `syntax/task/bash.py` parses bash scripts into `Assignment`/`Command` structures with interpolation analysis, but no pipeline stage uses it. The DAG stores the raw bash body as `script` with no static analysis.

**Proposed solution:** Wire `bash.Parser.parse()` into `Checker._load_import()` (for `.sh` files) and `force.py:_tool_definition()`. After parsing, validate:
- All variable references in interpolation (`${var}`) match known workflow inputs or task outputs.
- Shell syntax errors are caught early.
- The parsed `Script` is stored in the DAG's task definition alongside the raw body, so an executor can use it for fast variable lookup without re-parsing.

---

## Priority 3: Execution and serialization contract gaps

### P3.1: DAG JSON has no contract for evaluating interpolation defaults

**Spec:** String interpolation in defaults (`${outbase}.bam`) must be evaluated at execution time.

**Code:** `force.py:_build_task_param()` stores parsed defaults as `Word` → `Literal`/`Var`/`Expr` AST fragments serialized via `_interp_to_dict`. The executor receives `{"kind": "var", "name": "outbase"}` but there is no documented contract for how to resolve `Var` references against runtime bindings or evaluate `Expr` fragments. Examples and test coverage for the executor side don't exist.

**Proposed solution:** Define and document the interpolation JSON schema:
```json
{
  "kind": "word",
  "parts": [
    {"kind": "literal", "text": "path/to/"},
    {"kind": "var", "name": "sample_id"},
    {"kind": "literal", "text": ".bam"}
  ]
}
```
Also handle the case where a default is a single `Expr` (e.g., `${sample_id + "_out"}`) rather than just `Var` or `Literal`. Add a reference executor implementation that substitutes `Var` from a provided bindings dict and documents that `Expr` evaluation is deferred (the executor may choose to shell-evaluate or reject it). Document this in a `dag-schema.md` or in the DAG module docstring.

---

### P3.2: `pipe.swl` and `function.swl` produce different DAGs despite semantic equivalence

**Spec:** `A | B | C` desugars to a specific lambda form. Hand-written `Apply` + `Update` composition should produce the same DAG.

**Code:** `_normalize_output_value()` (`force.py:793`) handles the chain case's flat merge correctly (producing `Record` with all outputs as top-level fields). The equivalent hand-written form reaches forcing as a deeply nested `Merge` tree, and `_normalize_output_value` doesn't fully flatten it — it preserves nested `Merge` nodes under the `result` key.

**Proposed solution:** Fix `_normalize_output_value()` to recursively flatten `Merge` structures into a single `Record` whenever possible:
1. When both sides of a `Merge` are `Record`, merge their fields into one `Record`.
2. When one side is `Record` and the other is `Merge`, recursively normalize the `Merge` side and then merge fields.
3. Add `force.py:_canonicalize_merges()` called from `_force_apply` and `_final_outputs` to eagerly normalize `Merge` → `Record` at construction time, preventing the divergence from propagating.

---

## Priority 4: Parser correctness

### P4.1: Parser doesn't allow lambda on the right side of `|`

**Spec:** `chain ::= expr ( "|" expr)+` where `expr ::= ... | lambda`. Precedence: `->` > `|`. So `a | \x -> body` should parse as `a | (\x -> body)`.

**Code:** `_parse_chain_expr()` calls `_parse_update_expr()` for the right operand, which reaches `_parse_term()` — and `_parse_term` doesn't handle `\`. Lambda on the right of `|` is a parse error.

**Proposed solution:** In `_parse_chain_expr()`, change the right operand call from `_parse_update_expr()` to `_parse_expr()`. The precedence is already correct: `_parse_expr` handles `\` before delegating to `_parse_simple_expr`, and `->` binds tighter than `|` per spec, so `\x -> a | b` still parses as `\x -> (a | b)`.

---

### P4.2: Parser doesn't enforce spec's block grammar

**Spec:** `block ::= (binding eol)* expr` — only bindings before the final expression.

**Code:** `_parse_block()` allows any expression in any position, only validating that the final expression is not a `Binding`. `a = 1\nb\nc = 2\nd` parses without error despite `b` being a bare identifier in the middle.

**Proposed solution:** In `_parse_block()`, add a check after parsing each expression (except the last) that it is a `Binding`. Raise `ValueError` if a non-binding expression appears mid-block. This makes the parser match the spec grammar exactly.

---

## Priority 5: Code quality

### P5.1: `syntax/wf/node.py`: `Expr.__repr__` accesses subclass-only attributes

`__repr__` dispatches on `self.type` and accesses attributes (`self.name`, `self.body`, `self.fun`, etc.) that only exist on specific subclasses. Adding a new `NodeType` without updating `__repr__` causes `AttributeError` at display time.

**Proposed solution:** Replace the single `__repr__` with a per-subclass `__repr__` or a dispatcher dict that maps `NodeType` → format function. Use `getattr(self, field, '<missing>')` as a fallback for defensive access.

---

### P5.2: `ir/node.py`: `Unknown` sentinel can propagate to forcing

`lower.py` returns `ir.Unknown()` as a fallback for unhandled AST node types. `force.py:force_value()` has no handler for `ir.Unknown` and raises a generic `ValueError` without pointing to the original AST location or node type.

**Proposed solution:**
1. Replace the silent `return ir.Unknown()` fallback in `lower.py` with a `raise ValueError(...)` that names the unhandled `NodeType` and the expression's string representation.
2. Alternatively, add an `ir.Unknown` handler in `force_value()` that raises a descriptive error with the IR node's context.

---

### P5.3: `ir/dag.py`: `StepCall` and `MappedStep` serialization type confusion

`to_dict()` uses `getattr(step, 'map', None)` on all steps, relying on the attribute being absent on `StepCall`. `from_dict()` checks for the `map` key to decide which class to instantiate. If a `StepCall` somehow acquires map data, it's silently dropped during serialization.

**Proposed solution:** Use a single step class with an optional `map` field, or add an explicit `type` discriminator to `to_dict()` output that `from_dict()` uses regardless of key presence. If keeping two classes, use `@dataclass` inheritance with a discriminator field.

---

### P5.4: `ir/dag.py`: `_run_value_from_dict` called twice per item

The dict comprehension in `from_dict` calls `_run_value_from_dict` in both the value expression and the `if` filter, doubling execution. The function returns `None` for multiple semantically different reasons (param not found, value equals default, unparseable value) without distinguishing them.

**Proposed solution:** Replace the dict comprehension with a `for` loop that calls the function once, captures the result, and conditionally includes it. Return an explicit sentinel (e.g., `_SKIP = object()`) instead of `None` for the "skip this param" case, so `None` is reserved for actual errors.

---

### P5.5: `semantic/task/type.py`: `Param.__init__` type annotations mismatch

`Param.__init__` declares `typ: TypeKind`, `default: str`, `desc: str` (all non-optional), but `None` is regularly passed for all three throughout the codebase. This produces LSP type errors in every file that constructs `Param` with missing metadata.

**Proposed solution:** Change the annotation to `typ: TypeKind | None = None`, `default: str | None = None`, `desc: str | None = None`. Fix all call sites that already pass `None` explicitly.

---

### P5.6: `compile.py`: All exceptions handled identically

The `__main__` block catches all exceptions and prints the traceback to stdout with `os.EX_DATAERR`. There is no distinction between user errors (bad workflow input) and internal bugs. Traceback goes to stdout instead of stderr.

**Proposed solution:**
1. Print traceback to `sys.stderr` instead of stdout.
2. Define custom exception classes (`UserError` vs `InternalError`). User errors (bad syntax, type mismatch) print only the message; internal errors (unexpected `None`, attribute errors) print the full traceback.
3. Add `--verbose` flag to force full tracebacks for user errors during debugging.

---

### P5.7: `WorkflowCheck.chain_errors` and `WorkflowCheck.issues` are vestigial

Both properties in `WorkflowCheck` return `self.errors`. They are unused aliases from an earlier design that distinguished error categories.

**Proposed solution:** Remove both properties. Callers should access `check.errors` directly. If a future error-category distinction is needed, it should use a dedicated error type with a `category` field, not method aliases.
