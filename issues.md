# Issues

## ✅ Resolved

### Output type inference for Record outputs (Resolved)

`_flatten_nested_values` now recursively flattens all Records with dotted-key prefixes (e.g. `{a: {b: 3}}` → `{a.b: 3}`). Record values are eliminated from outputs before `_infer_output_type` is called, so the `Record → None` path is never reached.

### Output descriptions (Resolved)

`_output_desc` now propagates `desc` from `Input` and `Field(Input)` outputs, not just `Field(StepCall)`.

### Complex binding type inference (Resolved)

Covered by recursive flattening — only leaf types reach `_infer_output_type`.

### Transpiler table cleanup (Resolved)

WDL `map_by` is supported via `collect_by_key`. NF/SMK record-binding rejections documented as safety checks.

---

## 🟡 Open

### Lexer `+` token not in parser (Open)

The `+` token is now recognized by the lexer but the parser does not yet handle it. The spec does not define `+` semantics. Update the parser when the grammar is extended.

### Optionality inference for generated-map inputs

Core paths cover optionality for direct `_input()` calls. Generated-map inputs (created via `_emit_mapped_step`) may miss the `?` suffix from the source task's input types. The `input_schema` construction in `_emit_mapped_step` needs to preserve optionality.

## Known remaining gaps

| Gap | Status |
|-----|--------|
| Output descriptions (`OutputSpec.desc`) | **Fixed** — `Input`/`Field(Input)` desc now propagated via `_output_desc` |
| Output type inference for complex bindings (nested fields, records) | **Fixed** — recursive flattening eliminates Records before type inference |
| Optionality inference for generated-map inputs | Open — `_emit_mapped_step` `input_schema` may miss `?` suffix |
| WDL transpiler: no `map_by` (grouped scatter) support | **Resolved** — supported via `collect_by_key` |
| NF transpiler: rejects record bindings | Safety check; records flattened before transpilation |
| `+` operator token in lexer, not in parser | Needed if grammar extends to string concatenation / arithmetic |
