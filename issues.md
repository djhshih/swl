# Issues

## Output type inference for Record outputs (Open)

`_infer_output_type` returns `None` for `Record` values in workflow outputs (force.py). Transpilers default to `string`/`String`. Per spec, there is no composite type inference for records yet — this is a documented limitation rather than a bug.

## Lexer `+` token not in parser (Open)

The `+` token is now recognized by the lexer but the parser does not yet handle it. The spec does not define `+` semantics. Update the parser when the grammar is extended.

---

## Resolved

- **Output descriptions populated** — `_build_output_specs` now propagates `desc` from step `task.outputs[name].desc` for `Field(StepCall, name)` outputs.
- **Output type inference for nested Field chains** — `_infer_output_type` now walks nested `Field` chains to resolve the terminal step output's type.
- **NF transpiler record-binding error messages** — Improved to include field names and actionable guidance.
- **Optionality inference paths** — `_prune_unused_inputs` now routes through `_input()` method which checks `endswith('?')`.
- **WDL transpiler `map_by` support** — Implemented via `collect_by_key()` in prior session.
