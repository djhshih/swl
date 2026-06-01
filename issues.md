# Issues

Issues are ordered by priority. **Priority 1** items are architectural — they stem from the three-pass pipeline design (checker → lowerer → forcer) where analysis results computed by one stage are discarded before the next. **Priority 2** items are missing compile-time checks that violate the "fail fast" principle. **Priority 3** items are execution/serialization contract gaps. **Priority 4** items are parser correctness bugs. **Priority 5** items are code quality.

---

## Priority 1: Architectural — three-pass pipeline has no constraint bridge

The pipeline evaluates the program three times with different value domains:

| Pass | Input | Domain | Output |
|------|-------|--------|--------|
| Checker | WF AST | `OpenRecord`, `FunctionValue`, `ClosureValue`, `ComputationValue`, ... | `WorkflowCheck` (signature + type) |
| Lowerer | WF AST | `ir.Name`, `ir.Input`, `ir.Ref`, `ir.Variable`, `ir.Apply`, `ir.Lambda` | IR tree |
| Forcer | IR tree | `dag.Input`, `dag.Record`, `dag.Field`, `dag.StepCall`, `ForcedFunction` | `DAG` |

The Checker's analysis is not embedded in the IR. The Forcer re-derives everything from scratch. Three architectural issues follow from this.

---

### P1.1: Partial application model drops provenance between checker and forcer

**Status: Fixed.**

**Root cause:** Checker and forcer each evaluate partial application independently with different value domains. The checker produces `ClosureValue(function, bound_fields)` where `bound_fields` is `ClosedRecord` of `UnknownValue` entries — it remembers that fields are bound but not what they are bound to. The forcer produces `ForcedFunction(node, bound_value, signature)` with concrete bound values. There is no validation that the forcer's concrete evaluation agrees with the checker's abstract prediction.

**Spec violation:** "Partial application returns a new function." The checker's `ClosureValue` loses the bound values, so later stages may produce different step call bindings than the checker assumed. The checker also eagerly creates `ComputationValue` (treating a function as applied) when the bound record is an `OpenRecord`, making inference approximate.

**Fix:** Threaded the checker's satisfied-input information through the IR and into the forcer:
1. Added `satisfied: Set[str] = field(default_factory=set)` to `ir.Apply` (`ir/node.py`).
2. Checker records satisfied input sets per Apply AST node in `Checker._record_apply_satisfied()` (`semantic/wf/check.py`). For `ComputationValue` (full application), satisfied = all function input names. For `ClosureValue` (partial), satisfied = only explicitly-bound fields (where `_is_explicitly_bound_value()` returns True).
3. Lowerer reads `checker._apply_satisfied[id(expr)]` and passes as `satisfied=` to `ir.Apply` (`ir/lower.py:127`).
4. Added `satisfied: Set[str]` to `dag.ForcedFunction` (`ir/dag.py`) so partial applications carry accumulated satisfied provenance across multiple apply steps.
5. Forcer's `_force_apply()` extracts `satisfied` from `ir.Apply` and passes to `_apply()`. The `_apply()` method accumulates `fn.satisfied | satisfied` and stores it on new `ForcedFunction` instances for partial applications or passes it to emit methods for full applications (`ir/force.py`).
6. Added `Forcer._validate_bindings(inputs, satisfied, function_name)` that raises `ValueError` when step call bindings differ from the checker's predicted satisfied input set (`ir/force.py`). Called from both `_emit_task_call()` and `_emit_workflow_call()`.

---

### P1.2: Checker-inferred input metadata is lost between checker and DAG output

**Status: Fixed.** Implemented the thorough solution: added `ir.Input(name, type, desc)` as a first-class IR node. The lowerer creates `ir.Input` for names that the checker identified as demanded workflow inputs. The forcer handles `ir.Input` nodes directly, creating `dag.Input` with the metadata from the checker's `TaskSignature`.

**Root cause:** The checker infers `inputs = {a: file, b: str}` and produces a `TaskSignature`. The lowerer receives this signature and attaches it to the root `ir.Lambda`, but the forcer's `_input()` method (`force.py:723`) only receives a `spec` parameter when the caller happens to pass one. The `_refine_input_metadata()` recovery pass (`force.py:556`) scans step-level task input specs for matching names — input names that are only used as lambda parameters or record fields never get type/desc metadata.

**Fix:** (a) Added `ir.Input(name, type, desc)` to `ir/node.py`. (b) The `Lowerer` stores the checker's `TaskSignature` as `self.signature`; when lowering a `NodeType.id` that is not in env or imports but matches a signature input name, it emits `ir.Input(name, type, desc)` instead of `ir.Name`. (c) `force_value()` in `force.py` handles `ir.Input` by creating `dag.Input(name, type, desc)` directly, carrying the checker-inferred metadata forward to the DAG.

---

### P1.3: IR has no way to distinguish "workflow input" from "undefined variable"

**Status: Fixed.** Added `ir.Input(name, type, desc)` IR node. The lowerer creates `ir.Input` for checker-demanded workflow inputs and raises `ValueError` for undefined identifiers. The forcer handles `ir.Input` directly; `ir.Name` only looks up `ForceEnv` and raises an error if not found.

**Root cause:** The IR had no `ir.Input` node. The lowerer produced `ir.Name(expr.name)` for any identifier not found in its environment — whether that identifier is a workflow input parameter or a typo. The forcer couldn't tell the difference. When it saw `ir.Name("align_inpt")` not in the `ForceEnv`, its only option was to create a new `Input`, silently accepting a typo as a workflow input.

**Fix:** (a) Added `ir.Input(name, type, desc)` to `ir/node.py`. (b) The `Lowerer` stores the checker's `TaskSignature` as `self.signature`; when lowering a `NodeType.id` not in env/imports, it checks `signature.inputs` — if the name is a demanded workflow input, it emits `ir.Input(name, type, desc)`; otherwise it raises `ValueError(f'Undefined variable: {name}')`. (c) In `force.py`, `ir.Input` creates `dag.Input` with metadata; `ir.Name` only resolves through `ForceEnv` — unbound names raise `ValueError(f'Undefined variable during forcing: {name}')`. No implicit `Input` creation from `ir.Name`.

---

### P1.4: Duplicate pattern-matching for imports and builtins across checker and lowerer

**Status: Fixed.**

**Root cause:** Both `semantic/wf/check.py` and `ir/lower.py` independently implement pattern-matching for `import("path")`, `map f xs`, and `map_by f key xs`. The checker needs it for type inference; the lowerer needs it for IR construction. Any change to the AST structure or built-in detection rules must be mirrored in both files.

**Fix:** Created `syntax/wf/builtins.py` with three shared functions:
- `match_import(expr)` — returns path string or `None`
- `match_map(expr)` — returns `(function_expr, arg_expr)` or `None`
- `match_map_by(expr)` — returns `(key_expr, function_expr, arg_expr)` or `None`

Both `semantic/wf/check.py` and `ir/lower.py` now import `builtins` from `swl.syntax.wf` and call `builtins.match_*()` instead of their own `self._match_*()` methods. The individual `_match_import`, `_match_map`, and `_match_map_by` methods were removed from both classes. The lowerer's `_match_builtin` method (which is specific to lowering logic) remains but now delegates to `builtins.match_map` and `builtins.match_map_by` internally.

---

## Priority 2: Missing compile-time checks

### P2.1: No DAG circularity check

**Status: Fixed.**

**Spec Sec 3:** "DAG graph will be constructed and checked for circularity."

**Code:** `DAG` stores `steps` with dependency lists (`deps`) but has no cycle detection. `force.py` can produce a DAG with circular step dependencies (e.g., through self-referential bindings or mutually-recursive imports not caught by the import circularity check).

**Fix:** Added `DAG.validate()` (`ir/dag.py`) that runs a DFS-based cycle detection on the step dependency graph. Also validates that all dependency references point to existing step IDs. Called from `Forcer._finalize_dag()` before returning the DAG.

---

### P2.2: Run parameters with missing values silently dropped

**Status: Fixed.**

**Spec:** Run parameters (`cpu`, `memory`, `time`, `image`) are required for execution. If they have no default, the workflow input must provide them.

**Code:** `_normalize_task_run()` (`force.py:459`) projects run param names from the bound record and silently omits any name where `_project_field` returns `_SENTINEL` (i.e., the value is not available). The DAG compiles successfully but the executor receives no value for the required run parameter.

**Fix:** Modified `_normalize_task_run()` to check each run parameter after projection. If a parameter's field is not found in the bound record and the parameter has no `parsed_default`, raises `ValueError` with the parameter name and task path. Parameters with defaults may still be omitted.

---

### P2.3: Type compatibility only checked for explicit `|` chain syntax

**Status: Fixed.**

**Spec Sec 2:** Both pipeline (`a | b`) and record union (`r1 // r2`) must be type-checked.

**Code:** `Checker._check_chains()` only walks `NodeType.chain` nodes. No type checking is performed for `Apply`-based composition or `Update` expressions. `_merge_update_values` merges records without checking type compatibility of same-named fields.

**Fix:**
1. `TypeChecker` instance is now stored on the `Checker` as `self._type_checker` so it's accessible during expression evaluation.
2. In `Checker._eval_expr` for `NodeType.apply`, when both the function and argument resolve to `FunctionValue`s with names registered in the type checker, calls `TypeChecker.check_chain(left_name, right_name)` to validate output-to-input type compatibility.
3. In `_merge_arg_values`, overlapping record field types are checked via `_value_type` — known-type conflicts are detected (analogous to the existing `_merge_tables` check).

---

### P2.4: Import resolution only works at the top-level block

**Status: Fixed.**

**Spec:** `import: str -> fun` should be usable anywhere an expression is valid.

**Code:** `_load_imports()` only scans the outer block for `name = import("path")` patterns. An `import` inside a lambda or nested block is parsed but never resolved. It fails at forcing time with an unresolved `Function` node.

**Fix:**
1. In `Checker._eval_expr` for `NodeType.apply`: when `builtins.match_import(expr)` returns a path, the import is resolved inline via `self._load_import(stem, full_path)` and a `FunctionValue` is returned directly. Base directory is derived from `self._loading[-1]`.
2. In `Lowerer.lower_expr` for the same pattern: the import is resolved via `self._lower_inline_import(path)`, which uses `checker._load_import()` to produce an `ir.Function` node.

---

### P2.5: No glob validation for task output defaults

**Status: Fixed.**

**Spec:** Output parameter defaults are string expressions evaluated at execution time. Glob patterns are common (e.g., `bam/*.bam`).

**Code:** Output defaults are parsed as interpolation words but never validated. A malformed glob pattern (e.g., `[invalid`) is silently accepted.

**Fix:** Added `Forcer._validate_output_default_glob()` (`ir/force.py`) that validates glob bracket balance (`[...]` matching) in output parameter defaults. Called from `_build_task_param()` when `kind == 'out'` and the parameter has a default. Added `_interp_word_to_text()` helper to reconstruct the original text from interpolation AST for validation.

---

### P2.6: `bash.py` two-stage validation not integrated

**Status: Fixed.**

**Spec (Compile-time Checks):** Catch as many issues as possible statically.

**Code:** `syntax/task/bash.py` parses bash scripts into `Assignment`/`Command` structures with interpolation analysis, but no pipeline stage uses it. The DAG stores the raw bash body as `script` with no static analysis.

**Fix:** Wired `bash.Parser.parse()` into the pipeline with variable reference validation:
1. Added `_validate_bash_variables(parsed_body, input_names, context_name)` in `semantic/wf/check.py` that walks the parsed bash script, tracking defined variables (starting with task input names), and reports errors for any `${var}` reference that doesn't match a known task input or a previously assigned bash variable.
2. Called from `Checker._load_import()` for `.sh` files: raises `ValueError` if any unresolved variable references are found at import time.
3. Called from `force.py:_tool_definition()`: also validates at forcing time for defense-in-depth.
4. Added `task` and `parsed_body` fields to the `Import` class.

---

## Priority 3: Execution and serialization contract gaps

### P3.1: DAG JSON has no contract for evaluating interpolation defaults

**Status: Fixed.**

**Spec:** String interpolation in defaults (`${outbase}.bam`) must be evaluated at execution time.

**Code:** `force.py:_build_task_param()` stores parsed defaults as `Word` → `Literal`/`Var`/`Expr` AST fragments serialized via `_interp_to_dict`. The executor receives `{"kind": "var", "name": "outbase"}` but there is no documented contract for how to resolve `Var` references against runtime bindings or evaluate `Expr` fragments. Examples and test coverage for the executor side don't exist.

**Fix:**
1. Added interpolation JSON schema documentation as a module-level docstring in `ir/dag.py`, covering all four kinds (`word`, `literal`, `var`, `expr`) with resolution rules for each.
2. Added `dag.resolve_default(default, bindings=None)` — a reference executor implementation that takes a serialized interpolation dict and a bindings dict, substitutes `Var` references from bindings, passes `Literal` through, and preserves `Expr` fragments as-is (deferred to the executor). `Word` parts are concatenated after resolution.

---

### P3.2: `pipe.swl` and `function.swl` produce different DAGs despite semantic equivalence

**Status: Fixed.**

**Spec:** `A | B | C` desugars to a specific lambda form. Hand-written `Apply` + `Update` composition should produce the same DAG.

**Code:** `_normalize_output_value()` (`force.py:793`) handles the chain case's flat merge correctly (producing `Record` with all outputs as top-level fields). The equivalent hand-written form reaches forcing as a deeply nested `Merge` tree, and `_normalize_output_value` doesn't fully flatten it — it preserves nested `Merge` nodes under the `result` key.

**Fix:** Added `force.py:_canonicalize_merges(value)` that recursively normalizes `Merge` trees:
- Both sides `Record` → merge fields into one `Record`
- One side `Record` + other side `Merge` → recursively normalize the `Merge`, then merge
- Also normalizes inside `Record` fields and `Field` sources

Called from `_force_apply()` (after each application) and `_final_outputs()` (at output time), ensuring pipe-desugared and hand-written forms converge to the same DAG structure.

---

## Priority 4: Parser correctness

### P4.1: Parser doesn't allow lambda on the right side of `|`

**Status: Fixed.** `_parse_chain_expr()` right operand changed from `_parse_update_expr()` to `_parse_expr()`. Now `a | \x -> body` parses as `Chain(a, Function(\x -> body))`.

**Spec:** `chain ::= expr ( "|" expr)+` where `expr ::= ... | lambda`. Precedence: `->` > `|`. So `a | \x -> body` should parse as `a | (\x -> body)`.

**Root cause:** `_parse_chain_expr()` called `_parse_update_expr()` for the right operand, which reaches `_parse_term()` — and `_parse_term` doesn't handle `\`.

**Fix:** Changed the right operand call from `_parse_update_expr()` to `_parse_expr()` in `parser.py:202`. The precedence is already correct: `_parse_expr` handles `\` before delegating to `_parse_simple_expr`, and `->` binds tighter than `|` per spec, so `\x -> a | b` still parses as `\x -> (a | b)`.

---

### P4.2: Parser doesn't enforce spec's block grammar

**Status: Fixed.**

**Spec:** `block ::= (binding eol)* expr` — only bindings before the final expression.

**Code:** `_parse_block()` allows any expression in any position, only validating that the final expression is not a `Binding`. `a = 1\nb\nc = 2\nd` parses without error despite `b` being a bare identifier in the middle.

**Fix:** Modified `_parse_block()` to check that every expression before the final one is a `Binding`. The check fires after consuming the `eol` separator — if execution continues past the `eol`, the expression is mid-block and must be a binding. Final expression detection (eof or `bend`) happens before the check, so trailing-comment blocks (`x = 1\nx\n# done`) parse correctly.

---

## Priority 5: Code quality

### P5.1: `syntax/wf/node.py`: `Expr.__repr__` accesses subclass-only attributes

**Status: Fixed.** Replaced `__repr__` with a `getattr`-based dispatcher dict. Each `NodeType` maps to a format function that uses `getattr(self, attr, '<missing>')` for defensive attribute access. Unknown/undefined nodes fall back to `<expr type=??>`.

---

### P5.2: `ir/node.py`: `Unknown` sentinel can propagate to forcing

### P5.2: `ir/node.py`: `Unknown` sentinel can propagate to forcing

**Status: Fixed.**

1. Replaced `return ir.Unknown()` fallback in `lower.py:lower_expr()` with `raise ValueError(f'Unhandled AST node type during lowering: {expr.type}, expression: {expr}')` — this provides the original AST location and node type instead of silently producing an Unknown node.
2. Added an `ir.Unknown` handler in `force.py:force_value()` that raises a descriptive error, ensuring any remaining Unknown nodes (e.g., from `Block.result` default) are caught with a clear message rather than falling through to the generic catch-all.

---

### P5.3: `ir/dag.py`: `StepCall` and `MappedStep` serialization type confusion

**Status: Fixed.** Unified `StepCall` and `MappedStep` into a single `StepCall` class with an optional `map` field. Removed the `MappedStep` alias entirely — all references replaced with `StepCall` directly. Fixed `cwl/emit.py:_canonical_binding()` to check `getattr(value.source, 'map', None)` instead of class name string comparison.

---

### P5.4: `ir/dag.py`: `_run_value_from_dict` called twice per item

**Status: Fixed.** Replaced the dict comprehension with a walrus operator pattern (`if (value := func(...)) is not None`) so the function is called once per item. The original double call was both wasteful and misleading — the function may have side effects or expensive computation in future changes.

---

### P5.5: `semantic/task/type.py`: `Param.__init__` type annotations mismatch

**Status: Fixed.** Changed annotations to `typ: TypeKind | None = None`, `default: str | None = None`, `desc: str | None = None` in `Param.__init__`.

---

### P5.6: `compile.py`: All exceptions handled identically

**Status: Fixed.**
1. Traceback printed to `sys.stderr` instead of stdout.
2. Added `UserError` exception class — user errors (bad syntax, type mismatch) print only the message; internal errors print full traceback.
3. Added `--verbose` flag to force full tracebacks for user errors during debugging.

---

### P5.7: `WorkflowCheck.chain_errors` and `WorkflowCheck.issues` are vestigial

**Status: Fixed.** Removed both properties. Callers access `check.errors` directly.
