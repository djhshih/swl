# Plan: Address SWL Issues — Progress

## ✅ Completed

### 1. Output type inference for Record values (`dag/finalize.py:148-149`) — **Done** (plan2.md)

Replaced `_flatten_named_values` with `_flatten_nested_values` which recursively flattens Records using dotted-key prefixes. Eliminates all Records from output values before `_infer_output_type` is called, so the `Record` → `None` path is never hit.

**Changes:** `dag/finalize.py`, `transpile/common.py` (reconstruction helpers), all 4 transpilers, tests.

### 3. Output descriptions always `None` (`dag/finalize.py:99-104`, `_output_desc`) — **Done**

Added `Input` and `Field(Input)` branches to `_output_desc` to propagate description metadata from input parameter annotations.

**Changes:** `dag/finalize.py:_output_desc`.

### 4. Complex binding type inference — **Done** (covered by item 1)

Recursive flattening eliminates all nested Records/Merges from output values. After flattening, only leaf types (Input, Literal, Field) reach `_infer_output_type`, all of which are already handled. No additional inference logic needed.

### 6. Transpiler limitations — **Done** (docs update)

WDL `map_by` is already supported via `collect_by_key`. NF record-binding rejections are explicitly documented as safety checks.

### 7. Issues table consistency — **Done** (docs update)

Updated `issues.md` table to reflect current state.

---

## 🟡 In progress / Remaining

### 2. Lexer `+` token not handled by parser (`syntax/wf/lexer.py:220-221` vs `syntax/wf/parser.py`)

**Blocked on:** The spec does not define `+` semantics. Needs grammar extension first.

### 5. Optionality inference for generated-map inputs

`_emit_mapped_step` constructs `input_schema` from `_forced_signature` but needs to preserve the `?` suffix from source task input types for mapped-step inputs.
