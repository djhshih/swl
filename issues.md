# Issues

## Output type inference for Record outputs (Open)

`_infer_output_type` returns `None` for `Record` values in workflow outputs (force.py). Transpilers default to `string`/`String`. Per spec, there is no composite type inference for records yet — this is a documented limitation rather than a bug.

## Lexer `+` token not in parser (Open)

The `+` token is now recognized by the lexer but the parser does not yet handle it. The spec does not define `+` semantics. Update the parser when the grammar is extended.

## Known remaining gaps

| Gap | Status |
|-----|--------|
| Output descriptions (`OutputSpec.desc` is always `None`) | Not implemented — force.py always sets `desc=None` |
| Output type inference for complex bindings (nested fields, records) | Falls back to `None` → transpilers default to `string` |
| Optionality inference for all Input creation paths | Core paths covered; generated-map inputs may miss optionality |
| WDL transpiler: no `map_by` (grouped scatter) support | Explicitly rejected |
| NF transpiler: rejects record bindings | Safety check; records should be flattened before reaching transpilers |
| `+` operator token not in lexer | Needed if grammar extends to string concatenation / arithmetic |
