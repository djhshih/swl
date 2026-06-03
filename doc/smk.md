# Snakemake Transpiler: Implementation Status

## Overview

The Snakemake transpiler (`python/swl/transpile/smk/`) converts compiled SWL DAGs into Snakemake `Snakefile`-format output. It follows the same architecture as the CWL, WDL, and Nextflow transpilers. The transpiler is ~507 lines of Python (`emit.py`).

**Status:** Implemented. Covers simple task rules, pipelines, named ports, resources, config-based inputs, SWL output defaults, sub-workflow output defaults, shell expression interpolation, and mapped sub-workflows. `map_by` has a stub implementation.

---

## Architecture

```
transpile/smk/
├── __init__.py     # exports transpile_dag_dict, transpile_dag_file
├── __main__.py     # python -m entry point
├── cli.py          # CLI wrapper (uses transpile/_cli.py)
└── emit.py         # ~530 lines, the transpiler
```

### Entry points

- `transpile_dag_file(path)` — read DAG JSON from disk, transpile
- `transpile_dag_dict(data, workflow_id, _top_level, wrap_map)` — main transpilation from parsed DAG data

### Key internal functions

| Function | Purpose |
|----------|---------|
| `_task_to_rule(step, dag, wrap_map)` | Convert a task StepCall to a `rule {name}:` block with `input:`, `output:`, `params:`, `shell:`, and resource directives |
| `_binding_to_path(binding, port_name, is_mapped, scatter_ports, wrap_map)` | Convert a DAG binding value to a Snakemake expression for the `input:` section. Returns `config["name"]` for workflow inputs, SWL default spec patterns for inter-step references, and repr'd strings for Literal bindings. Only called for file-typed inputs |
| `_interpolate_shell(body, step)` | Rewrite `$var` / `${var}` / `${expr}` bash references to `{input.var}`, `{output.var}`, or `{params.var}`. Expression with operators (e.g. `${memory / cpu}`) becomes `$(( {params.memory} / {params.cpu} ))` |
| `_interp_to_smk(value)` | Render a SWL `word` interpolation dict into a Snakemake format string using `word_interp()` from `common.py` |
| `_collect_params(step)` | Gather Literal bindings, string-typed Input bindings, and `run` section values into `params:` entries |
| `_emit_resources(step)` | Map SWL `run` (cpu, memory, time, image) to Snakemake `threads:`, `resources:`, `container:` |
| `_dag_to_smk(dag, ...)` | Generate `rule all:` with target files for top-level workflow outputs |
| `_output_to_path(name, value, dag)` | Convert a workflow output value to a Snakemake path expression for `rule all:`. Uses `config["name"]` for Input values and SWL defaults for StepCall outputs |
| `_mapped_output_expand(name, step)` | Generate `expand("pattern", var=VAR)` for mapped workflow outputs, using the inner task step's output default pattern |
| `_collect_wildcard_globs(dag)` | Generate sweep variable assignments like `OUTBASE = config["outbase"]` from wildcards in mapped step output patterns |
| `_subworkflow_to_smk(step, parent_id)` | Recursively transpile sub-workflow steps (inline rules with `wrap_map` for mapped sub-workflows) |
| `_validate_supported(dag)` | Reject unsupported constructs: literal workflow outputs, Merge bindings |
| `_mapped_by_step_to_smk(step, dag, wrap_map)` | Stub for `map_by` grouped scatter — generates a two-rule structure with TODO grouping logic |

---

## SWL → Snakemake Mapping

### Well Supported

| SWL Concept | Snakemake Output | Details |
|---|---|---|
| **Task step** | `rule {name}:` with `shell:` | Shell body is interpolated: `$var` → `{input.var}` / `{output.var}` / `{params.var}` |
| **Pipeline chain (desugared)** | Sequential rules with matching output patterns | Dependencies are explicit via `Field(source=StepCall)` in bindings, rendered as SWL default output patterns |
| **Named input ports** | `input:` entries for files, `params:` for strings | File-typed inputs use `config["name"]` in `input:`; string/int/float inputs use `config["name"]` in `params:` |
| **Named output ports** | `output:` entries named by port | Output paths use the SWL task's `out` block `default` pattern (e.g., `{outbase}.bam`) |
| **`run.cpu`** | `threads: <value>` | Also added to `params:` for shell expression resolution |
| **`run.memory`** | `resources: mem_mb=<value>` | Also added to `params:` for shell expression resolution |
| **`run.time`** | `resources: runtime_minutes=<value>` | Also added to `params:` for shell expression resolution |
| **`run.image`** | `container: "docker://<value>"` | |
| **String interpolation (`$var`, `${var}`, `${expr}`)** | `{input.var}`, `{output.var}`, `{params.var}`, `$(( ... ))` | Simple vars: classified by name (input port → `{input.}`, output port → `{output.}`, param → `{params.}`). Expressions with operators like `${memory / cpu}` → `$(( {params.memory} / {params.cpu} ))` |
| **`word` interpolation in output defaults** | Snakemake format string `{var}` | `word_interp()` from `common.py` with `var_fn=lambda n: '{' + n + '}'` |
| **Workflow inputs** | `config["name"]` in `input:` section | User provides paths via Snakemake `--config` or `configfile:` |
| **Sub-workflow (no map)** | Inlined rules | Recursively transpiled; rules from sub-workflow appear in parent output |
| **Mapped sub-workflow (scatter)** | Inlined rules + sweep variables from config | `OUTBASE = config["outbase"]` at top of file; task outputs use SWL defaults; `rule all:` uses `expand()` |

### Weakly Supported

| SWL Concept | Status | Issue |
|---|---|---|
| **Mapped sub-workflow inter-step refs** | Works but fragile | Inter-step references use the prev step's SWL output default pattern (e.g., `{outbase}.bam`). If two steps have the same output pattern, Snakemake may see ambiguous rule matching |
| **Top-level map (not sub-workflow)** | Fallback path | `_default_output_path()` generates `results/{step_id}/{scatter_ports...}/{out_name}` — works but produces unwieldy paths with all scatter ports in directory structure |
| **`map_by` (grouped scatter)** | Stub only | Generates two rules (grouping + processing) but grouping logic is `# TODO` — no actual grouping implementation |
| **`cpu`/`memory` from `run` in shell** | Added as `{params.xxx}` | Snakemake convention is `{threads}` and `{resources.mem_mb}`. Works but non-idiomatic |
| **Expression-based output defaults** | `${{{...}}}` escape | SWL `${expr}` interpolations use double-brace Python escaping. Feasible but the transpiler hasn't been tested with expression defaults |

### Not Supported (explicitly rejected)

| SWL Concept | Behavior |
|---|---|
| **Literal workflow outputs** | Raises `ValueError`: "Snakemake does not support literal workflow outputs" |
| **Merge bindings** | Raises `ValueError`: "Merge bindings must be flattened before Snakemake transpilation" |
| **Record bindings** | Raises `ValueError`: "Record binding with fields [...] must be flattened before Snakemake transpilation" |
| **Record workflow outputs** | Flattens into individual target files via recursive `_output_to_path()` call — experimentally supported but not tested |

### Snakemake-specific details

- **Config-based inputs**: File-typed inputs use `config["port_name"]` in the `input:` section. String/int/float-typed inputs use `config["port_name"]` in the `params:` section. Shell references use `{input.port}` for files and `{params.port}` for strings.
- **Output path resolution**: Task rule outputs use the SWL task's `out` block `default` pattern rendered through `_interp_to_smk()`. For mapped sub-workflows, sweep variables (`OUTBASE = config["outbase"]`) are generated from wildcards found in the output patterns by traversing the sub-workflow DAG to find inner task steps' output defaults. Non-mapped sub-workflows also traverse the inner DAG via `_inner_output_default()`.
- **Shell interpolation**: Regex `(?<!\$)\$\{(.+?)\}|(?<!\$)\$(\w+)` matches `${var}`, `${expr}`, and `$var` but not `$${escaped}`. Simple vars are replaced with `{input.var}`, `{output.var}`, or `{params.var}`. Expressions with operators (e.g. `${memory / cpu}`) become `$(( {params.memory} / {params.cpu} ))` using bash arithmetic.
- **Scatter via config lists**: For mapped sub-workflows, the sweep variable is defined as `VAR = config["var"]` and the `rule all:` target uses `expand("pattern", var=VAR)`. The user provides a list of identifiers in config (e.g., `config["outbase"] = ["sample1", "sample2"]`).

---

## Support Matrix

| Feature | Status | Tests |
|---------|--------|-------|
| Simple task rule | ✅ Well supported | function, pipe, explicit |
| Pipeline chain (|) | ✅ Well supported | function, pipe (equivalence) |
| Named I/O ports | ✅ Well supported | All |
| Resources (cpu/memory/time/image) | ✅ Well supported | All |
| Config-based inputs (file→`input:`, string→`params:`) | ✅ Well supported | All |
| SWL output default patterns | ✅ Well supported | function, pipe, explicit |
| String interpolation (`$var`, `${var}`, `${expr}`) | ✅ Well supported | function, pipe, explicit |
| Sub-workflow output defaults (mapped & non-mapped) | ✅ Well supported | partial, import_partial, map, panel |
| Mapped sub-workflow (scatter) | ⚠️ Weak support | map, panel |
| Top-level map step | ⚠️ Weak support | panel (merge input) |
| `map_by` (grouped scatter) | ❌ Stub only | map_by |
| Literal workflow outputs | ❌ Rejected | — |
| Merge bindings | ❌ Rejected | — |
| Record bindings | ❌ Rejected | — |

### Coverage

The test suite (`bash test.sh`):
- Generates `tests/smk/*.smk` from all compiled DAGs
- Compares `pipe.smk` and `explicit.smk` against `function.smk` for equivalence
- `map.smk`, `panel.smk`, `map_by.smk` are generated but have no golden comparison (they exercise the weaker support paths)

---

## Known Limitations

1. **Inter-step Snakemake ambiguity**: When two task steps use the same output default pattern (e.g., both `align` and `sort` have `bam='{outbase}.bam'`), Snakemake may be unable to determine which rule produces a given output file. The SWL model relies on explicit port wiring; Snakemake relies on pattern uniqueness. Users should ensure distinct output patterns across steps.

2. **`map_by` grouping logic**: Only the rule structure is emitted. The actual grouping logic (reading column data, partitioning by key, writing per-group files) is a `# TODO` placeholder.

3. **Top-level mapped steps**: The `_default_output_path()` fallback includes ALL scatter ports in the directory path, producing unwieldy paths like `results/{step}/{port1}/{port2}/.../{portN}/{output}`.

