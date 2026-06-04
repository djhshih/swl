# Plan: Address SWL Issues

## 1. Output type inference for Record values (`dag/finalize.py:148-149`)

**Problem:** `_infer_output_type` returns `None` for `Record` values. Then `DAG._check_output_types` raises `ValueError` because `OutputSpec.type` is `None`.

**Root cause:** `_flatten_named_values` in `finalize.py:67-75` only does one level of Record expansion via `_flatten_merge_value`. A nested Record (e.g. `{b: {c: 1, d: 2}}`) survives into `_build_output_specs` as a Record sub-value.

**Fix plan:**
1. In `_flatten_named_values`, recurse until no more Record values remain (exhaustive flattening) — or handle in `_flatten_outputs`.
2. As a safety net, in `_infer_output_type`, add a `Record` branch that infers the type from field values. If all fields share the same type, use it; otherwise fall back to `"str"`.
3. Update tests to cover nested-record and mixed-field Record outputs.

---

## 2. Lexer `+` token not handled by parser (`syntax/wf/lexer.py:220-221` vs `syntax/wf/parser.py`)

**Problem:** The lexer emits `TokenType.plus` (line 220-221) but the parser ignores it — no parse method handles `+`, no AST node exists for it.

**Blocked on:** The spec does not define `+` semantics. Needs grammar extension first.

**Fix plan (when spec is updated):**
1. Add `NodeType.plus` to `syntax/wf/node.py` and a `Plus` node class.
2. Add `_parse_plus_expr` in `parser.py` at an appropriate precedence level (likely between `update` and `apply`).
3. Add binary-expression handling in the evaluator/checker/lowerer.

---

## 3. Output descriptions always `None` (`dag/finalize.py:87-92`, `_output_desc`)

**Problem:** `_output_desc` only handles `Field(StepCall).` For `Input` and `Record` outputs the description is lost even when metadata exists.

**Fix plan:**
1. Add `isinstance(normalized, Input)` branch: propagate `normalized.desc`.
2. Add `isinstance(normalized, Field)` + `isinstance(normalized.source, Input)` branch: propagate `normalized.source.desc`.
3. Add `isinstance(normalized, Record)` branch: try unwrapping single-field Records to peek at the inner value's desc. For multi-field Records, descriptions are not meaningful at the whole-record level.
4. **Note:** `tooldefs.py:15` already sets `desc` on Input from task param specs; `_refine_input_metadata` (`finalize.py:56-64`) merges input descs. These already work — the gap is purely in `_output_desc`.

---

## 4. Complex binding type inference (nested fields, Records)

**Problem:** `_infer_output_type` returns `None` for `Record` and has limited handling for complex `Merge`/`Field` chains.

**Fix plan (same as item 1 fix, extended):**
1. Record inference: infer type from field values (same as fix 1.2).
2. Merge inference: `_normalize_output_value` collapses record-record merges into a single Record; handle the resulting Record (falls into fix 1.2). Non-record merges are rejected during DAG validation anyway, so no inference needed.
3. Field chains (e.g. `Field(Field(StepCall, "a"), "b")`): the existing `_root_field_source` + `_step_field_type` path in `_infer_output_type` (lines 155-159) already handles 2+ deep field chains.

---

## 5. Optionality inference for generated-map inputs

**Problem:** `_emit_mapped_step` (`dag/emit.py:63-84`) constructs `input_schema` from `_forced_signature` but does not explicitly carry optionality (`?`) from the source task's input specs into mapped-step inputs.

**Fix plan:**
1. In `_emit_mapped_step`, when building `input_schema`, preserve the `?` suffix from `normalize_swl_type(spec.type.value)` for each input that has an optional type in the task signature.
2. Trace callers of `_emit_mapped_step` to ensure the optionality flag is set on mapped input ports.
3. Verify that `_input()` in `tooldefs.py:12-18` and `_refine_input_metadata` already handle `?`-ending types — if so, the fix may only need to ensure the type string from the task signature (with `?`) flows through to the mapped step's `input_schema`.

---

## 6. Transpiler limitations (documented)

**WDL `map_by`:** Already supported via `collect_by_key` (`transpile/wdl/emit.py:357-418`). The issues table is outdated — mark as resolved.

**NF record bindings:** Explicitly rejected (`transpile/nf/emit.py:272-277`). Records must be flattened before reaching the NF transpiler. This is a safety check, not a bug. Document clearly.

---

## 7. Issues table consistency

**Problem:** The "Known remaining gaps" table in `issues.md` has an outdated row: "`+` operator token not in lexer". The `+` token IS now in the lexer (it's the parser that doesn't handle it). Also, the WDL `map_by` row is outdated.

**Fix:** Update the table to match current reality:
- Change "`+` operator token not in lexer" → "`+` operator token in lexer, not in parser".
- Mark WDL `map_by` row as resolved.
- Mark NF record bindings row as "Intentionally rejected (safety check)".
