# Implementation Plan

## Phase 1: Workflow parser - COMPLETE

The workflow parser is already structured well as three modules:
- `python/swl/syntax/wf/lexer.py`
- `python/swl/syntax/wf/node.py`
- `python/swl/syntax/wf/parser.py`

That structure should be the model for task parsing as well, with one important difference:
- workflow syntax needs a lexer because it is expression-oriented
- task annotation syntax is primarily line-oriented

## Phase 2: Task parser refactor

### Current state

Task parsing is currently split across:
- `python/swl/syntax/task/annotation.py`
- `python/swl/syntax/task/parser.py`
- `python/swl/syntax/task/type.py`

But the implementation is incomplete and inconsistent:
- `annotation.py` mixes parsing logic with ad-hoc data structures
- `parser.py` is an older incomplete parser sketch and currently contains syntax/runtime issues
- `type.py` contains semantic/type-checking concepts, not syntax nodes
- there is no clean separation between:
  - extracting annotation text from shell comments
  - parsing annotation syntax
  - parsing interpolation syntax
  - representing parsed task metadata as structured nodes
  - semantic typing/checking

### Target structure

Refactor task handling into separate syntax and semantic layers.

Recommended syntax-side files:
- `python/swl/syntax/task/node.py`
- `python/swl/syntax/task/parser.py`
- `python/swl/syntax/task/interpolation.py`

Recommended semantic-side file:
- `python/swl/semantic/task/type.py`

This separates:
- task annotation syntax
- interpolation syntax
- semantic typing/checking

## Phase 2.1: Define task AST/data nodes in `task/node.py`

`task/node.py` should represent parsed task syntax, not type-checker logic.

Recommended nodes:

- `Task`
  - `annotation: Annotation`
  - `body: str`

- `Annotation`
  - `doc: str | None`
  - `sections: list[Section]`
  - or, if simpler for downstream use:
  - `inputs: list[Param]`
  - `outputs: list[Param]`
  - `run: list[Param]`

- `Section`
  - `kind: SectionType`
  - `params: list[Param]`

- `Param`
  - `names: list[str]`
  - `type: str | None`
  - `default: Interpolation | None`
  - `desc: str | None`

Suggested enum:
- `SectionType = Enum('SectionType', ['in_', 'out', 'run'])`

Important design choice:
- keep these nodes declarative
- do not mix type compatibility or workflow inference logic into `task/node.py`
- interpolation values can be stored as parsed interpolation objects from `task/interpolation.py`

## Phase 2.2: Where should `type.py` live?

You renamed `node.py` to `type.py`, which is an improvement because that file is not a node module.

Recommendation:
- yes, move it to `python/swl/semantic/task/type.py`

Reason:
- the contents are semantic/type-system logic, not syntax
- `TypeKind`, `TaskSignature`, `TypeChecker`, and compatibility checks belong to semantic analysis
- keeping them under `syntax/` blurs the boundary between parsing and semantic validation

So the better long-term ownership is:
- `python/swl/syntax/task/node.py` -> parsed task syntax structures
- `python/swl/syntax/task/parser.py` -> task annotation parser
- `python/swl/syntax/task/interpolation.py` -> interpolation parser + interpolation nodes
- `python/swl/semantic/task/type.py` -> type definitions and type checking

If you want a minimal migration path:
1. leave `python/swl/syntax/task/type.py` temporarily in place
2. add a new `python/swl/semantic/task/type.py`
3. move imports gradually
4. eventually delete or shim the old location

## Phase 2.3: Implement `task/parser.py` as the annotation parser

`annotation.py` should be split so that:
- parsing lives in `task/parser.py`
- parsed data structures live in `task/node.py`

Recommended responsibility split:

- `task/parser.py`
  - public entry point: `Parser.parse(script_content: str) -> node.Task`
  - split shell script into annotation region and body region
  - normalize comment prefixes (`#` and optional following space)
  - parse task doc line (`@ ...`)
  - parse sections: `in`, `out`, `run`
  - parse parameter lines, including:
    - multiple names: `a, b file`
    - optional type: `file?`
    - defaults: `= 2`, `= "x"`, `= ${outbase}.bam`
    - descriptions: `| text`
    - continuation description lines beginning with `|`
  - call `interpolation.py` to parse default values where needed
  - return both parsed annotation and raw body

- `task/node.py`
  - data classes / node classes only

Then `annotation.py` can either:
- be removed entirely, or
- become a thin compatibility wrapper during transition

## Phase 2.4: Parsing strategy for annotation syntax

Unlike workflow syntax, task annotations are line-oriented.

A simple modular parser can work in two layers:

### Layer 1: comment extraction / normalization

Input is full shell script text.

Produce:
- normalized annotation lines
- raw body text

By:
- reading top-of-file comment lines that form the annotation block
- removing leading `#`
- removing one optional following space
- preserving line order
- preserving blank lines where meaningful
- stopping annotation capture when shell code begins

This layer is separate from syntax parsing.

### Layer 2: line-based annotation parser

Consume normalized annotation lines with a small cursor-based parser.

Suggested helpers in `task/parser.py`:
- `_at(i=0)` -> current normalized line
- `_eat()` -> consume line
- `_eof()` -> no more lines
- `_parse_doc()`
- `_parse_section()`
- `_parse_param()`
- `_parse_names()`
- `_parse_type()`
- `_parse_default()`
- `_parse_desc()`

This keeps the structure similar in spirit to `wf/parser.py`, even though it operates on lines rather than tokens.

## Phase 2.5: Do we need a lexer like in `wf/`?

For annotation parsing: probably not at first.

Why a lexer is not strictly needed for annotations:
- annotation syntax is mostly line-oriented, not expression-oriented
- section boundaries are line-based (`in`, `out`, `run`)
- parameter syntax is shallow enough to parse directly from each line
- defaults can be handed off to `interpolation.py` rather than fully tokenized by the annotation parser

Recommendation for annotations:
- do not build a full lexer yet
- first implement a clean line-based `task/parser.py`
- if parameter parsing becomes messy, introduce a very small param-line lexer later

So for task annotations the likely end state is:
- definitely `task/node.py`
- definitely `task/parser.py`
- probably no separate annotation lexer yet

## Phase 2.6: Bash interpolation should be a separate shared parser module

The task file contains at least two syntaxes:
- annotation comments at the top
- interpolation syntax used in defaults and shell body text

These should not be handled by one parser.

Reason:
- annotation parsing is metadata parsing
- interpolation parsing is shell-like expression parsing
- combining them will make the task parser harder to reason about and harder to test

Examples from `tests/*.sh`:
- `${outbase}.bam`
- `${outbase}.bai`
- `${cpu}`
- `${memory / cpu}`

These appear in both:
- annotation defaults, e.g. `#   bam file = ${outbase}.bam`
- shell body commands, e.g. `> ${outbase}.bam`

Recommendation:
- keep interpolation in one shared file for now:
  - `python/swl/syntax/task/interpolation.py`

Responsibilities of `interpolation.py`:
- define interpolation nodes
- parse interpolation fragments/words
- parse variable references: `$x`, `${x}`
- parse concatenated interpolated words like `${outbase}.bam`
- conservatively represent braced expressions like `${memory / cpu}`
- be reusable from both:
  - `task/parser.py` when parsing annotation defaults
  - future bash-body analysis when scanning shell code

Important:
- this does not need to be a full bash parser initially
- it should parse only the subset needed by current task files and tests
- interpolation logic should be implemented once here and reused everywhere else

## Phase 2.7: Keep `interpolation.py` as one file for now

To keep the refactor manageable, `interpolation.py` should contain both:
- interpolation node definitions
- interpolation parsing logic

That means no separate `interpolation_node.py` for now.

Suggested shapes in `python/swl/syntax/task/interpolation.py`:
- `Word(parts)`
- `Literal(text)`
- `Var(name)`
- `Expr(text)`

Examples:
- `${outbase}.bam` -> `Word([Var('outbase'), Literal('.bam')])`
- `${cpu}` -> `Word([Var('cpu')])`
- `${memory / cpu}` -> `Word([Expr('memory / cpu')])`

This is a good compromise:
- modular enough to keep interpolation separate from annotation parsing
- simple enough to avoid too many tiny files right now

## Phase 2.8: Suggested parser layering for task files

Recommended flow:

1. `task/parser.py`
   - split file into annotation region and body region
   - parse annotation region into `Annotation`
   - keep body as raw text
   - parse annotation defaults through `interpolation.py`
   - return `Task(annotation, body)`

2. `task/interpolation.py`
   - parse interpolation values from annotation defaults
   - also serve as the shared interpolation parser for shell body analysis
   - later, a bash-body analyzer can call into this module rather than reimplement interpolation parsing

3. optional future `task/bash.py`
   - scan/analyze shell body
   - reuse `interpolation.py` for `${...}` and `$x` parsing

4. `semantic/task/type.py`
   - build task signatures
   - perform type validation
   - handle compatibility checks for workflows/tasks

This keeps concerns separated and allows tests for each layer independently.

## Phase 2.9: Concrete recommendation on parser boundaries

Yes, interpolation should be a different parser/module from the annotation parser.

We should have:
- one parser for task annotations
- one shared parser for interpolation syntax
- one semantic module for task typing/checking

Not necessarily a full bash parser, but definitely a separate interpolation parser from the annotation parser.

The key rule is:
- parse interpolation once in `task/interpolation.py`
- reuse it from both annotation parsing and later bash-body analysis

## Phase 2.10: Proposed concrete APIs

### Annotation/task parser API

- `Parser.parse(script: str) -> node.Task`
- optionally `parse_file(path: str) -> node.Task`

Behavior:
- returns parsed task annotation plus raw body
- raises `ValueError` on malformed annotation syntax
- parses annotation defaults into interpolation nodes

### Interpolation parser API

In `python/swl/syntax/task/interpolation.py`:
- `parse_word(s: str) -> Word`
- maybe `parse_braced(s: str) -> Var | Expr`

Behavior:
- parses interpolation fragments used in defaults and shell words
- raises `ValueError` on malformed interpolation syntax

### Semantic type API

In `python/swl/semantic/task/type.py`:
- `parse_type(type_str: str) -> TypeKind`
- `types_compatible(output_type, input_type) -> bool`
- `TaskSignature`
- `TypeChecker`

## Phase 2.11: Validation boundaries

Validation to include in `task/parser.py`:
- task doc line must begin with `@`
- at least one section must exist
- only `in`, `out`, `run` are valid section names
- parameter lines must contain at least one name
- malformed annotation defaults should raise `ValueError`

Validation to include in `interpolation.py`:
- malformed `${...` expressions
- malformed variable references

Validation to defer to semantic/type phase:
- duplicate parameter names across sections
- compatibility with workflow chaining
- required/optional semantics beyond surface syntax
- interpretation of arithmetic-like shell expressions inside `${...}`

## Phase 2.12: Suggested migration steps

1. Keep `wf/` untouched
2. Add `python/swl/syntax/task/node.py` for actual task syntax/data nodes
3. Rewrite `python/swl/syntax/task/parser.py` from scratch as a clean line-based parser
4. Add `python/swl/syntax/task/interpolation.py` as a single file containing interpolation nodes + parser
5. Move current `python/swl/syntax/task/type.py` to `python/swl/semantic/task/type.py`
6. Update imports gradually or leave a temporary compatibility shim
7. Turn `task/annotation.py` into either:
   - a compatibility shim, or
   - remove it after imports are updated
8. Add focused unit tests for annotation parsing and interpolation parsing

## Phase 2.13: Unit tests needed for annotation parsing

Create parser tests analogous to `wf/parser.py` tests.

At minimum cover:
- task doc only + one section
- `in`, `out`, `run` sections
- multiple params in one section
- multiple names on one line
- optional types (`file?`)
- array types (`[file]`)
- defaults with bare words
- defaults with quoted strings
- defaults with interpolation `${outbase}.bam`
- inline descriptions with `|`
- continuation description lines
- malformed section header
- malformed parameter line
- missing task doc or missing section, depending on desired strictness

## Phase 2.14: Unit tests needed for interpolation parsing

Add tests for:
- `${outbase}.bam`
- `${outbase}.bai`
- `${cpu}`
- `${ref}`
- `${memory / cpu}`
- plain literals with no interpolation
- mixed literal + variable parts
- malformed `${...` expressions
- quoted strings if we decide to support them immediately

## Phase 2.15: Recommended file ownership after refactor

### Syntax layer
- `python/swl/syntax/wf/lexer.py`
  - workflow tokenization only
- `python/swl/syntax/wf/node.py`
  - workflow AST only
- `python/swl/syntax/wf/parser.py`
  - workflow parser only

- `python/swl/syntax/task/node.py`
  - task syntax/data nodes only
- `python/swl/syntax/task/parser.py`
  - task annotation parser only
- `python/swl/syntax/task/interpolation.py`
  - interpolation nodes + interpolation parser
- `python/swl/syntax/task/annotation.py`
  - temporary compatibility wrapper or remove

### Semantic layer
- `python/swl/semantic/task/type.py`
  - `TypeKind`, `TaskSignature`, `TypeChecker`, compatibility helpers

## Immediate recommendation

Implement the task parser in a modular way without introducing a lexer for annotations.

Specifically:
- add `task/node.py` for parsed task syntax/data
- rewrite `task/parser.py` as a line-oriented annotation parser
- keep `task/interpolation.py` as a single shared file for interpolation nodes + parser
- use `task/interpolation.py` from annotation defaults now and from bash-body analysis later
- move type-checking logic to `swl/semantic/task/type.py`
- treat annotation parsing, interpolation parsing, bash-body analysis, and semantic typing as separate concerns
