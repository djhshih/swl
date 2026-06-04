# Shell Variable Interpolation

## Overview

SWL task scripts use bash-style `${var}` substitutions that resolve to each target workflow language's native expression syntax. The interpolation infrastructure is shared across all four transpilers (SMK, WDL, NF, CWL).

## Shared Infrastructure

### Parsing

Bash variable references are parsed by `python/swl/syntax/task/interpolation.py`:

- `parse_word(s)` — converts `${var}`, `$var`, `${expr}`, and literal text into a `Word([parts...])` tree
- Node types: `Literal(text)`, `Var(name)`, `Expr(text)`

The bash parser (`python/swl/syntax/task/bash.py`) uses `parse_word` for each word in a bash script body, producing a `Script` AST of `Command` and `Assignment` statements.

### Shared utilities in `python/swl/transpile/common.py`

**`word_interp(value, literal_fn, var_fn, expr_fn, joiner)`** — renders a single annotation-style interpolation dict (from DAG JSON) by applying per-part callbacks. Used for output path defaults.

**`interp_script(body, var_fn, expr_fn)`** — renders a raw bash script body by regex-matching `${var}`, `${expr}`, and `$var` and applying per-match callbacks. Returns the fully interpolated script string. Used for command-block interpolation.

**`classify_var(name, input_names, output_names, run_names)`** — determines a variable's scope:
| Scope | Meaning | Resolution |
|---|---|---|
| `'input'` | Task input port | Native input syntax |
| `'run'` | Task run param (cpu, memory, time) | Native run-param syntax or inline literal |
| `'output'` | Task output port | Native output syntax |
| `'shell'` | Bash builtin (`$HOME`, `$PATH`, etc.) | Pass through verbatim |
| `None` | Unknown / not in scope | Pass through verbatim |

### Bash built-in allowlist

Defined in `python/swl/syntax/task/bash.py` as `_BUILTIN_VARS`:

```python
_BUILTIN_VARS = frozenset({
    'HOME', 'PATH', 'USER', 'SHELL', 'PWD', 'RANDOM',
    'LINENO', 'SECONDS', 'UID', 'HOSTNAME', 'OSTYPE',
    'BASH', 'BASH_VERSION', 'BASH_ENV', 'IFS', 'PS1',
    'PS2', 'PS3', 'PS4', 'OLDPWD', 'SHLVL', 'TERM',
    'LANG', 'LC_ALL', 'LC_CTYPE', 'DISPLAY', 'TMPDIR',
    'EDITOR', 'VISUAL', 'PAGER',
})
```

Builtins are skipped by the variable checker and passed through verbatim by all transpilers.

### Variable checker (`python/swl/semantic/wf/bashvars.py`)

`_validate_bash_variables(parsed_body, known_vars, context)` walks bash Script statements sequentially, maintaining a mutable `defined` set (initialized from `known_vars` = `input_names | run_names`). For each:
- `Assignment`: checks RHS var refs against `defined`, adds LHS to `defined`
- `Command`: checks all var refs in each word against `defined`

Both `Var` and `Expr` parts are checked. Builtins are excluded from validation.

## Per-Transpiler Behavior

### SMK (`python/swl/transpile/smk/emit.py`)

| Pattern | Resolution |
|---|---|
| `${var}` / `$var` (input port) | `{input.var}` |
| `${var}` / `$var` (output port) | `{output.var}` |
| `${var}` / `$var` (param/run) | `{params.var}` |
| `${expr}` | `$(( resolved_expr ))` — inner variables resolved per scope |
| `$HOME`, `$PATH`, etc. | Pass through verbatim |
| Unknown var | Pass through verbatim |

**Scope**: Input names exclude string/int/float types (these go in `params:`). Output names are from `step.outputs`. Run params are `cpu`, `memory`, `time`.

**Implementation**: `_interpolate_shell(body, step)` → `interp_script(body, var_fn, expr_fn)` with scope-specific `classify_var`-style lookup.

### WDL (`python/swl/transpile/wdl/emit.py`)

| Pattern | Resolution |
|---|---|
| `${var}` (known var) | `~{var}` |
| `${expr}` (known vars inside) | `~{expr}` |
| `$var` (unbraced) | `{input.cpu}` / `~{cpu}` / `${task.cpus}` / `inputs.cpu` |
| `$HOME`, `$PATH`, etc. | Pass through verbatim |
| Unknown var | Pass through verbatim (`${var}` unchanged) |

**Scope**: `known_vars = set(input_names) | run_var_names`. Run vars (`cpu`, `memory`, `time`) with values are added as WDL task inputs with defaults, so `~{cpu}` and `~{memory / cpu}` resolve correctly in WDL.

Both `$var` (unbraced) and `${var}` (braced) are handled identically — they are equivalent bash syntax.

**Implementation**: `_interpolate_bash_vars(body, known_vars)` → calls `interp_script(body, var_fn, expr_fn)`.

### NF (`python/swl/transpile/nf/emit.py`)

| Pattern | Resolution |
|---|---|
| `${var}` (input port) | `${var}` (pass through — Nextflow template resolves from `input:` scope) |
| `${cpu}` | `${task.cpus}` (mapped via `_NF_RUN_MAP`) |
| `${memory}` | `${task.memory}` |
| `${time}` | `${task.time}` |
| `${expr}` with run vars | `${task.cpus}`, `${task.memory}` inside expression |
| `$HOME`, `$PATH`, etc. | Pass through verbatim |
| Unknown var | Pass through verbatim |

**Scope**: Input names from `task.inputs`. Run names mapped through `_NF_RUN_MAP = {'cpu': 'cpus', 'memory': 'memory', 'time': 'time'}` to `${task.*}` syntax.

**Implementation**: `_interpolate_shell(body, step)` → `interp_script(body, var_fn, expr_fn)`.

### CWL (`python/swl/transpile/cwl/emit.py`)

| Pattern | Resolution |
|---|---|
| `${var}` (file input) | `inputs.var.path` (in CWL `$(...)` JS expression) |
| `${var}` (string input) | `inputs.var` |
| `${cpu}` (run param) | Inlined literal value (e.g., `2`) |
| `${memory}` (run param) | Inlined literal value (e.g., `8192`) |
| `${expr}` | `(resolved_expr)` — inner vars replaced per scope |
| `$HOME`, `$PATH`, etc. | Pass through verbatim as JS string literals |
| Unknown var | Pass through verbatim |

**How it works**: The script body is converted to a CWL `$(...)` JavaScript expression that constructs the shell script string at runtime. The `InitialWorkDirRequirement` `entry` becomes:
```json
{
  "entryname": "script.sh",
  "entry": "$('bwa mem -t ' + 2 + ' ' + inputs.ref.path + ' > ' + inputs.outbase + '.bam')"
}
```

String inputs → `inputs.name`, file inputs → `inputs.name.path`, run params → inlined JavaScript values. `InlineJavascriptRequirement` is automatically added to the tool requirements when the body contains variable references.

**Multi-line scripts**: Each line becomes a JS expression joined with `"\\n"`:
```json
"entry": "$('samtools sort -@ ' + 2 + ' aligned.bam ' + inputs.outbase + '.bam' + '\" + \"\\n\" + \"' + 'samtools index ' + inputs.outbase + '.bam ' + inputs.outbase + '.bai')"
```

**Implementation**: `_interpolate_shell(body, step)` builds the JS expression directly (does not use `interp_script`), using `classify_var()` from `common.py`.

## Resolution Summary

| ${var} scope | SMK | WDL | NF | CWL |
|---|---|---|---|---|
| Input (file) | `{input.ref}` | `~{ref}` | `${ref}` | `inputs.ref.path` |
| Input (string) | `{params.outbase}` | `~{outbase}` | `${outbase}` | `inputs.outbase` |
| Run (cpu) | `{params.cpu}` | `~{cpu}` | `${task.cpus}` | Inlined `2` |
| Run (memory) | `{params.memory}` | `~{memory}` | `${task.memory}` | Inlined `8192` |
| Run (time) | `{params.time}` | `~{time}` | `${task.time}` | Inlined `30` |
| Output | `{output.bam}` | N/A* | `${bam}` | N/A* |
| Shell builtin | Verbatim `$HOME` | Verbatim `$HOME` | Verbatim `$HOME` | Verbatim `"$HOME"` |
| Unknown | Verbatim `${unknown}` | Verbatim `${unknown}` | Verbatim `${unknown}` | Verbatim `"${unknown}"` |

\* Output ports are not valid references in command blocks at runtime.
