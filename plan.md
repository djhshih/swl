# Remaining implementation plan

## Goal

Drive `python/swl/` toward production readiness. Tracks concrete implementation work derived from:
- `dag.md` — the DAG JSON spec and cross-transpiler requirements
- `cwl.md` — CWL transpiler status and gaps
- `wdl.md` — WDL 1.1 transpiler implementation plan
- `nf.md` — Nextflow transpiler implementation plan
- `spec.md` — SWL language specification
- `progress.md` — current gap list

---

## Priority 5: Remaining quality issues

Driven by `issues.md`. Lower-urgency items that affect correctness or completeness but don't block basic compilation.

### 5a. Output type inference for complex bindings

**Status:** DONE. `_infer_output_type` in `force.py` now walks nested `Field` chains to resolve the terminal step output's type. `Record` outputs still return `None` (per spec — no composite type inference yet).

Tests: `test_infer_output_type_nested_field_chain_resolves_terminal_step`, `test_infer_output_type_nested_field_chain_mapped_resolves_to_array`, `test_infer_output_type_record_returns_none`.

### 5b. Output descriptions on workflow outputs

**Status:** DONE. `_build_output_specs` now propagates `desc` from step `task.outputs[name].desc` when the output value is `Field(StepCall, name)`, `Input`, and `Literal` sources stay `None`.

Tests: `test_output_desc_propagated_from_step_output`, `test_output_desc_none_for_input_source`, `test_output_desc_none_for_literal_source`.

### 5c. Optionality inference for generated-map inputs

**Status:** DONE. `_prune_unused_inputs` now routes map-by source input creation through `_input()` method (which applies the `?` suffix check), instead of constructing `Input(...)` directly.

Tests: `test_map_by_source_input_created_through_input_method`.

### 5d. WDL transpiler: `map_by` (grouped scatter) support

**Status:** DONE (completed in prior session). `_mapped_by_step_to_wdl` emits `collect_by_key()` calls with proper `Pair` types. See `tests/dag/map_by.json` and `tests/wdl/map_by.wdl`.

### 5e. NF transpiler: graceful record-binding rejection

**Status:** DONE. Error messages now include field names and guidance: `"Record binding with fields [...] must be flattened before Nextflow transpilation. ... Use explicit field bindings instead."`. Added compile-time guard in `_finalize_dag`.

Tests: `test_record_binding_error_includes_fields`.

### 5f. Lexer: `+` operator token

**Status:** DONE. `plus` added to `TokenType` enum and recognized in lexer as a single-character token.

---

## Dependency graph

```
Q5a (output type inference) ──── done
Q5b (output descriptions) ────── done
Q5c (optionality gaps) ───────── done
Q5d (WDL map_by) ─────────────── done
Q5e (NF record error) ────────── done
Q5f (lexer +) ────────────────── done
```

---

## File change summary

| File | What changes |
|------|-------------|
| `ir/force.py` | (Q5a) extend `_infer_output_type` for nested Field chains; (Q5b) propagate step output `desc` to `OutputSpec.desc`; (Q5c) route map-by source input creation through `_input()` |
| `ir/test_force.py` | Tests for Q5a, Q5b, Q5c, Q5e (3 + 3 + 1 + 1 = 8 new tests) |
| `transpile/wdl/emit.py` | (Q5d) implement grouped-scatter (map_by) support — completed in prior session |
| `transpile/wdl/test_emit.py` | Tests for Q5d — completed in prior session |
| `transpile/nf/emit.py` | (Q5e) improve error messages for record bindings (both `Record` instance and `source='record'` dict paths) |
| `transpile/nf/test_emit.py` | Tests for Q5e |
| `syntax/wf/lexer.py` | (Q5f) recognize `+` token — added `plus` to TokenType enum and single-character handling |
