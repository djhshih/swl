# Shell Variable Interpolation Fix Plan

## Problem

Each SWL task's shell script contains `${ref}`, `${cpu}`, `${outbase}.bam` etc. These are bash variable substitutions that need to be resolved to each target workflow language's native expression syntax. Currently:

| Transpiler | Resolution | Status |
|---|---|---|
| **SMK** | `{input.ref}` / `{params.cpu}` / `{output.outbase}` | ✅ Correct |
| **WDL** | `~{ref}` / `~{cpu}` / `~{outbase}` | ❌ `cpu` not in WDL input scope |
| **NF** | verbatim `${ref}` / `${cpu}` | ❌ `${cpu}` is undefined in Nextflow |
| **CWL** | verbatim `${ref}` / `${cpu}` | ❌ undefined bash variables at runtime |

## Current Infrastructure

### Shared interpolation nodes (`python/swl/syntax/task/interpolation.py`)
- `Word([parts...])`, `Literal(text)`, `Var(name)`, `Expr(text)` — the core node types for both annotation and bash parsing
- `parse_word(s)` converts `${var}`, `$var`, `${expr}`, and literal text into these nodes
- Used by: annotation parser (`parser.py:110`), bash parser (`bash.py:64,70`), and the bash variable checker (`bashvars.py:6`)

### Annotation interpolation flow

```
SWL source:  # out bam file = ${outbase}.bam
                         ↓
parser.py:110  interpolation.parse_word("${outbase}.bam")
                         ↓
               Word([Var('outbase'), Literal('.bam')])
                         ↓
tooldefs.py:213  _interp_to_dict(param.default)
                         ↓
               {"kind": "word", "parts": [
                 {"kind": "var", "name": "outbase"},
                 {"kind": "literal", "text": ".bam"}]}
                         ↓
               Stored in DAG JSON as step.task.outputs.<name>.default
                         ↓
               Each transpiler reads the JSON dict and applies
               word_interp() with transpiler-specific callbacks:
                 SMK: _interp_to_smk(value)   →  {outbase}.bam
                 WDL: _interp_to_wdl(value)   →  ~{outbase}.bam
                 NF:  _interp_to_nf(value)    →  ${outbase}.bam
                 CWL: _interp_to_cwl_glob(v)  →  $(inputs.outbase + '.bam')
```

Key: `word_interp()` in `common.py:44` is the shared generic function — it takes a JSON dict + 4 callback functions and returns the rendered string. All four transpilers use it for annotation defaults.

### Bash body interpolation flow

```
SWL source:  bwa mem -t ${cpu} ${ref}
                         ↓
bash_parser.py:66  _parse_words() → per-word: interpolation.parse_word(part)
                         ↓
               [Word([Var('cpu')]), Word([Var('ref')])]
                         ↓  (stored in Script.statements[i].words)
               Used by bashvars.py for validation ONLY.
                         ↓
               NOT used by transpilers! Instead each ad-hoc re-parses:
                 SMK: regex  r'(?<!\$)\$\{(.+?)\}|(?<!\$)\$(\w+)'
                 WDL: brace-matching state machine
                 NF:  nothing (verbatim)
                 CWL: nothing (verbatim)
```

**Key problem**: The bash parser already produces the correct parse tree (Var/Expr/Literal), but transpilers ignore it and re-parse with ad-hoc methods. The annotation flow shares this code correctly; the bash body flow does not.

### Variable checker (`python/swl/semantic/wf/bashvars.py`)
- `_validate_bash_variables(parsed_body, known_vars, context)` walks statements sequentially maintaining a mutable `defined` set
- `known_vars` = `set(input_names) | set(run_names)` from task annotation
- `defined` starts as a copy of `known_vars`, then grows as assignments are encountered
- For each `Assignment` statement:
  1. Check all `Var` refs in RHS against `defined` (error if unregistered)
  2. Add LHS name to `defined` (now available for subsequent statements)
- For each `Command` statement:
  1. Check all `Var` refs in each word against `defined` (error if unregistered)

Example:
```bash
tmp=aligned.bam                # Assignment: checks RHS (refs: none), adds 'tmp' to defined
samtools sort ${tmp} ${outbase} # Command: checks 'tmp' ✓, 'outbase' ✓ (from input_names)
echo ${nonexistent}             # Command: 'nonexistent' NOT in defined → ERROR
```

**Limitation**: Only checks `Var` parts, not `Expr` parts (e.g. `${memory / cpu}`). The `memory` and `cpu` refs inside the expression text are never validated.

### Checker invocation
- `tooldefs.py:138-142`: `known_vars = set(input_names) | set(run_names)` → validate
- `imports.py:46-49`: same pattern
- Both sites call `_validate_bash_variables(parsed_body, known_vars, ...)` at DAG construction time, before transpilation

### Interpolation per transpiler

- **SMK** `_interpolate_shell`: regex-based, resolves `${var}`/`$var` by looking up in input/output/param sets
- **WDL** `_interpolate_bash_vars`: brace-matching state machine, converts `${var}` → `~{var}` with NO scope resolution
- **NF**: no interpolation — script body passes verbatim
- **CWL**: no interpolation — script body passes verbatim in InitialWorkDirRequirement

## Streamlining approach

The annotation interpolation and bash body interpolation share the same underlying parser (`interpolation.parse_word`) and should share the same rendering approach. Currently annotation interpolation reuses `word_interp()` in all four transpilers, but bash body interpolation uses ad-hoc parsers (or nothing).

**Fix**: Add `_interpolate_script(script, word_interp_fn)` to `transpile/common.py` that walks a parsed `Script` and applies `word_interp_fn` to each word in each command. This mirrors how `word_interp` works for single words but operates on the full script AST.

Each transpiler would then:
1. Call `bash.parse(body)` (once, shared)
2. Define per-word callbacks for `word_interp` that handle scope resolution
3. Call `_interpolate_script(parsed_script, per_word_callbacks)` to produce the rendered body
4. Remove the ad-hoc regex/state-machine parsers

This gives us:
- **Consistent parsing**: All transpilers use the same `interpolation.parse_word()`, correctly distinguishing `Var('cpu')` from `Expr('memory / cpu')`
- **Removes code duplication**: No more regex in SMK, no more state machine in WDL
- **Annotation + bash parity**: Both use `word_interp`-based callbacks; only the word extraction differs (annotation: single word from DAG JSON; bash: many words from parsed Script)

## Plan improvements

### Issue 1: No bash built-in allowlist

The checker (`bashvars.py`) would reject `$HOME`, `$PATH`, `$RANDOM` or other well-known bash variables because they're not in the `defined` set. These are valid bash references that should pass through verbatim, not be resolved.

**Fix**: Add a `_BUILTIN_VARS` frozen set to `bash.py` or `bashvars.py`:

```python
_BUILTIN_VARS = frozenset({
    'HOME', 'PATH', 'USER', 'SHELL', 'PWD', 'RANDOM',
    'LINENO', 'SECONDS', 'UID', 'HOSTNAME', 'OSTYPE',
    'BASH', 'BASH_VERSION', 'IFS', 'PS1',
})
```

The checker skips validation for builtins. Transpilers leave `${HOME}` verbatim.

### Issue 2: Parsed Script re-parsed per transpiler (wasteful)

Currently `step.task['body']` is a raw string. `tooldefs.py:138` already constructs `parsed_body = task_bash.Parser().parse(task.body)` during DAG construction but discards it after validation — it's never stored.

Each transpiler would call `bash.parse(body)` independently → 4 redundant parses. Moreover, the DAG `step.task['body']` reaches transpilers as a raw string in the JSON dict.

**Fix**: Either:
- (A) Store `parsed_body` on the `StepCall` object during DAG construction, so transpilers receive it alongside `step.task`
- (B) Add `'parsed_body'` to the DAG JSON dict (as a serializable format)
- (C) Simplest: let transpilers call `bash.parse()` once each — the cost is negligible for typical script sizes

Option (C) is pragmatic. But for Phase 3's `interp_script` to work, we need to agree on the interface: does each transpiler pass a `Script` object or a raw string? If raw string, `interp_script` parses internally and every transpiler benefits without changing the DAG pipeline.

### Issue 3: No transpile-time safety check for unresolved vars

After interpolation, there's no check for leftover `${...}` patterns in the output. If the transpiler's resolver misses a variable (e.g., a typo or an unexpected scope), the generated code silently contains an unresolved reference.

**Fix**: Add `_check_no_remaining_vars(body)` to each transpiler's shell interpolator that scans for `${...}` patterns after resolution and raises an error if any remain. This catches bugs in the transpiler itself (not user errors — those are caught by the checker at DAG construction).

### Issue 4: CWL approach is underspecified

The plan says "use EnvVarRequirement or inline assignments" but both options have problems:

- `EnvVarRequirement` entries take literal values or CWL expressions. Run params (`cpu`, `memory`) are **not** CWL expression-accessible — they live in `ResourceRequirement`, not `inputs`. So you can't write `$(resources.cpu)` in an env var definition.
- Inline assignments (`cpu=2`) work for literal-known values (run params, literal bindings) but not for dynamic inputs (workflow inputs, step output snapshots).

CWL's fundamental constraint: `InitialWorkDirRequirement` entries are flat files. CWL has no template engine for file content. Variables can only be injected via:
1. `EnvVarRequirement` — but only supports literals and `$(inputs.*)` expressions
2. Converting the `entry` to a CWL expression that constructs the script string at runtime
3. Adding extra input ports to each tool for run-param values, then referencing them as `inputs.var`

**Recommended fix**: Approach (2). Convert the InitialWorkDir `entry` from a literal string to a CWL `$(...)` expression. The expression uses `inputs.var` for input variables and inlines literal values for run params. Memory values (already in MB) and cpu values (already int) are inlined directly. This gives correct CWL v1.0:

```json
{
  "class": "InitialWorkDirRequirement",
  "listing": [{
    "entryname": "script.sh",
    "entry": "$('samtools sort -@ ' + 2 + ' -m ' + (8192 / 2) + ' aligned.bam ' + inputs.outbase + '.bam')"
  }]
}
```

Pros: Valid CWL, no new input ports, no env var injection. Cons: `file`-typed variable names (like `ref`) become `inputs.ref.path` but the original script expects just the filename. CWL stages files to cwd with their original basenames, so `inputs.ref.path` gives the full path — this may differ from the expected bare filename.

### Issue 5: No interpolation-specific golden file tests

The plan mentions "update golden files" but doesn't specify adding test cases designed to catch interpolation regressions. Without them, a future change that accidentally reverts to verbatim pass-through would go undetected.

**Fix**: Add dedicated test DAGs/SWL files per transpiler that exercise known interpolation patterns:
- `${var}` simple variable
- `$var` unbraced (should same as `${var}`? document the decision)
- `${expr / with / vars}` expression
- `${memory / cpu}` resource expression
- Variables from each scope (input, output, run, assignment)
- Negative test: a var that should NOT be resolved (e.g., bash built-in `$HOME`)

Each produces a golden file that is compared in `test.sh`.

### Issue 6: SMK backward compatibility

Replacing SMK's regex `_interpolate_shell` with `interp_script()` may produce slightly different output even for the same inputs (e.g., whitespace handling, edge cases with `$` in the script). The existing golden file comparisons (`test.sh`) may fail.

**Fix**: When replacing SMK's interpolator, run `bash test.sh` immediately and diff golden outputs. If differences are cosmetic (whitespace), update golden files. If semantic (actual interpolation behavior changes), investigate.

### Issue 7: Phase ordering

The current ordering (Phase 1 → 2 → 3 → 4 → 5 → 6 → 7) intermixes transpiler-specific changes with shared infrastructure. Better to:

- Phase A: bash parser + checker + common.py helpers (shared infra)
- Phase B: All 4 transpiler changes together (consistent approach, cross-checkable)
- Phase C: Tests + golden file updates (after all behavior is settled)

### Issue 8: Scope resolver responsibility split

Currently the plan has each transpiler implementing its own `classify_var` + syntax mapping. This duplicates the classification logic (input vs run vs output vs shell) 4 times.

**Fix**: `common.py`'s `classify_var` returns the scope. Then each transpiler only needs a mapping table: `{scope: syntax_template}`. For example:

```python
# common.py
def classify_var(name, input_names, output_names, run_names):
    if name in input_names:
        return 'input'
    if name in run_names:
        return 'run'
    if name in _BUILTIN_VARS or name.startswith('_'):
        return 'shell'
    return 'unknown'

# per transpiler (e.g., NF)
_SYNTAX_MAP = {
    'input': '${name}',
    'run': '${task.nf_name}',
    'shell': '${name}',
    'unknown': None,  # leave verbatim
}
# with a per-transpiler name-remapping dict for run params: {cpu: cpus, memory: memory}
```

This centralizes the classification and makes each transpiler a thin syntax layer.

## Revised Plan

### Phase A: Shared infrastructure

1. **bash.py**: Add `_BUILTIN_VARS`, `iter_var_refs()`, `_extract_expr_vars()`
2. **bashvars.py**: Check `Expr` parts, skip builtins
3. **common.py**: Add `classify_var()` and `interp_script()`

### Phase B: All 4 transpiler changes

4. **SMK**: Replace regex with `interp_script()`
5. **WDL**: Replace brace-matcher with `interp_script()`; add run vars as optional task inputs
6. **NF**: Add `interp_script()` with scope resolution
7. **CWL**: Add `interp_script()` with CWL expression-based script construction

### Phase C: Tests and golden files

8. Add interpolation-specific test cases per transpiler
9. Generate golden files, update `test.sh`
10. Full test run — verify all transpilers produce correct output

## File-by-file change summary

| File | Change |
|------|--------|
| `python/swl/syntax/task/bash.py` | Add `_BUILTIN_VARS`, `iter_var_refs()`, `_extract_expr_vars()` |
| `python/swl/semantic/wf/bashvars.py` | Check `Expr` parts; skip builtins |
| `python/swl/transpile/common.py` | Add `classify_var()` and `interp_script()` |
| `python/swl/transpile/smk/emit.py` | Replace regex `_interpolate_shell` → `interp_script()` |
| `python/swl/transpile/wdl/emit.py` | Replace brace-matcher → `interp_script()`; add run vars as optional inputs |
| `python/swl/transpile/nf/emit.py` | Add `_interpolate_shell` via `interp_script()` (input → `${var}`, run → `${task.*}`) |
| `python/swl/transpile/cwl/emit.py` | Add `_interpolate_shell` via `interp_script()`; CWL expression `entry` for InitialWorkDir |
| `tests/unit/swl/syntax/task/test_bash.py` | Add `iter_var_refs` and builtins tests |
| `tests/unit/swl/transpile/wdl/test_emit.py` | Update `test_bash_var_interpolation` |
| `tests/unit/swl/transpile/nf/test_emit.py` | Add shell interpolation test |
| `tests/unit/swl/transpile/cwl/test_emit.py` | Add shell interpolation test |
