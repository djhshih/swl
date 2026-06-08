# Issues

## Open

### `+` operator token in lexer, not in parser

The `+` token is recognized by the lexer but the parser does not yet handle it. The spec does not define `+` semantics. Update the parser when the grammar is extended.
