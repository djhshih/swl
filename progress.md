# Progress

## Known remaining gaps

| Gap | Status |
|-----|--------|
| Output descriptions (`OutputSpec.desc` is always `None`) | Not implemented — force.py always sets `desc=None` |
| Output type inference for complex bindings (nested fields, records) | Falls back to `None` → transpilers default to `string` |
| Optionality inference for all Input creation paths | Core paths covered; generated-map inputs may miss optionality |
| WDL transpiler: no `map_by` (grouped scatter) support | Explicitly rejected |
| NF transpiler: rejects record bindings | Safety check; records should be flattened before reaching transpilers |
| `+` operator token not in lexer | Needed if grammar extends to string concatenation / arithmetic |
