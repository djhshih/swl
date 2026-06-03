# WDL Transpiler: Implementation Status

## Overview

The WDL transpiler (`python/swl/transpile/wdl/`) converts compiled SWL DAGs into WDL 1.1 `task` + `workflow` blocks. It follows the same architecture as the CWL, Nextflow, and Snakemake transpilers. The transpiler is ~520 lines of Python (`emit.py`).

**Status:** Implemented. Covers simple tasks, pipelines, named ports, requirements (cpu/memory/time/image), output path templates, shell variable interpolation (`${var}` → `~{var}`), scatter for mapped steps, `collect_by_key` for `map_by`, struct definitions for Record bindings, and inline sub-workflow blocks. `$var` (unbraced) interpolation and `${memory / cpu}` expression resolution are not fully handled.

---

## Architecture

```
transpile/wdl/
├── __init__.py     # exports transpile_dag_dict, transpile_dag_file
├── __main__.py     # python -m entry point
├── cli.py          # CLI wrapper (uses transpile/_cli.py)
└── emit.py         # ~520 lines, the transpiler
```

### Entry points

- `transpile_dag_file(path)` — read DAG JSON from disk, transpile
- `transpile_dag_dict(data, workflow_id, _top_level)` — main transpilation from parsed DAG data

### Key internal functions

| Function | Purpose |
|----------|---------|
| `_task_to_wdl(step)` | Convert a task StepCall to a WDL `task { }` block with `input {}`, `command <<< >>>`, `output {}`, `requirements {}` |
| `_interpolate_bash_vars(line)` | Convert `${var}` bash references to `~{var}` WDL interpolation in command blocks |
| `_interp_to_wdl(value)` | Render a SWL `word` interpolation dict into a WDL string template `"~{var}"` |
| `_emit_requirements(step)` | Map SWL `run` (cpu, memory, time, image) to WDL `requirements {}` |
| `_dag_to_wdl(dag, workflow_id, tasks)` | Generate the `workflow { }` block with `input {}`, `call` statements, and `output {}` |
| `_binding_to_wdl_expr(binding, current_step_id, dag)` | Convert a DAG binding to a WDL expression |
| `_mapped_step_to_wdl(step, tasks)` | Emit `scatter` block for mapped steps, handling per-row indexing vs broadcast inputs |
| `_mapped_by_step_to_wdl(step, tasks)` | Emit `collect_by_key()` + `scatter` block for `map_by` grouped scatter |
| `_subworkflow_to_wdl(step, parent_id)` | Recursively transpile sub-workflow steps (inline `workflow { }` blocks) |
| `_collect_structs(dag)` | Walk all Record bindings to generate `struct` definitions |
| `_emit_struct(record)` | Emit a WDL `struct` definition from a Record binding's field types |

---

## SWL → WDL 1.1 Mapping

### Well Supported

| SWL Concept | WDL 1.1 Output | Details |
|---|---|---|
| **Task step** | `task { }` with `command <<< >>>` | Script body embedded with `${var}` → `~{var}` interpolation |
| **Named input ports** | `input {}` block | `File`, `String`, `Int`, `Float`, `Array[File]` per SWL type |
| **Named output ports** | `output {}` with `File name = "~{var}.ext"` | Output paths use SWL `default` pattern rendered as `~{var}` template |
| **`run.cpu`** | `requirements { cpu: <value> }` | |
| **`run.memory`** | `requirements { memory: "<value> MB" }` | Value converted to MB |
| **`run.time`** | `requirements { time_minutes: <value> }` | |
| **`run.image`** | `requirements { container: "<value>" }` | |
| **Shell `${var}` interpolation** | `~{var}` in `command <<< >>>` | Brace-delimited bash vars (`${cpu}`) converted to WDL `~{cpu}`. WDL resolves from task input scope |
| **`word` interpolation in output defaults** | WDL string template `"~{var}"` | `word_interp()` from `common.py` with `var_fn=lambda n: f"~{{{n}}}"` |
| **Workflow inputs** | `workflow { input { ... } }` | Same type mapping as task inputs |
| **Input binding (workflow input)** | `<name>` | Direct reference to workflow input |
| **Input binding (literal)** | Inline literal value | Strings quoted, numbers bare, booleans `true`/`false` |
| **Step output binding** | `<alias>.<output_name>` | References task's named output |
| **Field projection** | `<expr>.<field_name>` | Member access on struct/pair |
| **`[file]` array types** | `Array[File]` | Properly mapped via `to_wdl_type()` |
| **`[str]`, `[int]`, `[float]` array types** | `Array[String]`, `Array[Int]`, `Array[Float]` | Properly mapped |
| **Sub-workflow (no map)** | Inline `workflow { }` block | Recursively transpiled; appears as separate workflow block, called from parent |
| **Record binding** | `struct` literal | `struct` definition emitted at top of file, binding emits `StructName { fields... }` |
| **Workflow output** | `output { Type name = <expr> }` | Type inferred from source binding; `Array[]` for mapped step outputs |

### Weakly Supported

| SWL Concept | Status | Issue |
|---|---|---|
| **`$var` (unbraced) in shell** | Not interpolated | The `_interpolate_bash_vars` function only converts `${var}` (with braces), not `$var` (without braces). `$cpu` passes through literally |
| **`${memory / cpu}` expression** | Passes through to `~{memory / cpu}` | The brace matcher converts `${memory / cpu}` → `~{memory / cpu}`, but `memory` and `cpu` are not WDL task input variables (they're in `requirements {}`). Would cause WDL evaluation error |
| **`run.cpu`/`memory` access in command** | Not accessible as `~{cpu}` | WDL requirements values (`cpu`, `memory`) are not automatically in scope for `command` blocks. `~{cpu}` would fail — needs `~{...}` to reference task inputs, not requirements |
| **Mapped sub-workflow scatter inputs** | All ports scatter-indexed | When a mapped sub-workflow has many scatter ports, all are indexed with `[i]` in the scatter body. Broadcast inputs are not distinguished from scatter inputs |
| **`map_by` deeply nested Pair types** | Correct but verbose | Tables with many columns produce deeply nested `Pair[..., Pair[...]]` types (e.g., 10 levels for 10 columns), making the WDL output hard to read |

### Not Supported (explicitly rejected)

| SWL Concept | Behavior |
|---|---|
| **Merge bindings** | Raises `ValueError`: "Merge bindings must be flattened before WDL transpilation" |
| **Literal workflow outputs** | Raises `ValueError`: "WDL does not support literal workflow outputs" |

### WDL-specific details

- **Shell interpolation**: The `_interpolate_bash_vars` function uses a brace-matching state machine to convert `${var}` → `~{var}`. It handles nested braces (e.g., `${foo_{$bar}}`). Only brace-delimited `${var}` is converted — `$var` (unbraced) passes through literally, which is fine for WDL because `$var` in `command <<< >>>` is treated as a shell variable, not WDL interpolation.
- **Output path templates**: SWL `word` interpolation defaults are rendered to WDL string expressions `"~{var}.ext"` using `_interp_to_wdl()`. The `File` declaration includes the template as the default value.
- **`requirements {}` vs `runtime {}`**: All resource specifications (`cpu`, `memory`, `time`, `image`) go into `requirements {}` (WDL 1.1 canonical). The `runtime {}` section is not emitted.
- **Sub-workflow call aliasing**: Sub-workflows are called with an `as` alias matching the step ID to avoid name collisions with the generated workflow name.
- **`map_by` source decomposition**: When `map.source` is `input`, the transpiler decomposes the source input into individual array columns declared in the workflow `input {}`. Column types are derived from `input_schema`.
- **Struct generation**: Record bindings are collected into WDL `struct` definitions at the top of the file. Field types are inferred from the binding source.

---

## Support Matrix

| Feature | Status | Tests |
|---------|--------|-------|
| Simple task | ✅ Well supported | function, pipe, explicit |
| Pipeline chain (|) | ✅ Well supported | function, pipe (equivalence) |
| Named I/O ports | ✅ Well supported | All |
| Requirements (cpu/memory/time/image) | ✅ Well supported | All |
| Workflow input/output | ✅ Well supported | All |
| Output path templates (`"~{var}"`) | ✅ Well supported | function, pipe, explicit |
| Shell `${var}` interpolation | ✅ Well supported | All |
| Literal input bindings | ✅ Well supported | partial |
| Field projections | ✅ Well supported | All |
| `[file]` / `[str]` / `[int]` / `[float]` types | ✅ Well supported | All |
| Sub-workflow (no map) | ✅ Supported | partial, import_partial |
| Record bindings / struct | ✅ Supported | — |
| Mapped step (scatter) | ⚠️ Weak support | map, panel |
| `map_by` (collect_by_key) | ⚠️ Weak support | map_by |
| `$var` (unbraced) interpolation | ❌ Not converted | — |
| `$memory / cpu` expression in command | ❌ Not in WDL scope | — |
| Merge bindings | ❌ Rejected | — |
| Literal workflow outputs | ❌ Rejected | — |

### Coverage

The test suite (`bash test.sh`):
- Generates `tests/wdl/*.wdl` from all compiled DAGs
- Compares `pipe.wdl` and `explicit.wdl` against `function.wdl` for equivalence of golden output
- `panel.wdl`, `map.wdl`, `map_by.wdl` are generated but have no golden comparison

---

## Known Limitations

1. **`$var` (unbraced) interpolation**: The `_interpolate_bash_vars` function only converts `${var}` (with braces) to `~{var}`. Unbraced `$var` references in command blocks pass through literally, which WDL's `command <<< >>>` passes to the shell as-is. This is correct for shell variables (`$PATH`, `$?`), but SWL variables like `$cpu` or `$ref` would not be interpolated by WDL.

2. **`~{memory / cpu}` expression**: `${memory / cpu}` is converted to `~{memory / cpu}` by the brace matcher. However, `memory` and `cpu` are defined in `requirements {}`, not as WDL task inputs. WDL `command <<< >>>` blocks only see task `input {}` variables in scope. The expression `~{memory / cpu}` would cause a WDL evaluation error. A fix would need to either add these as optional inputs or use WDL's `runtime` accessor if available.

3. **Mapped sub-workflow broadcast vs scatter**: When a mapped sub-workflow has scatter ports, all workflow inputs matching scatter port names are indexed with `[i]`. Inputs that should be broadcast (shared across all iterations) are incorrectly indexed as well. This only works correctly when all inputs are scatter ports.

4. **`map_by` deeply nested types**: Tables with many columns produce deeply nested `Pair[...]` types for `collect_by_key`. While functionally correct, the generated WDL is verbose and hard to debug. The `_col_access_path` function generates long `.right.right.left` accessor chains.

5. **Sub-workflow module isolation**: Sub-workflows are inlined as separate `workflow { }` blocks with duplicated input declarations. WDL 1.1 allows multiple workflow blocks (only the first is the entry point), but some engines may warn or behave unexpectedly with multi-workflow files.

---

## Example Output

### Simple pipeline (`function.swl`)

Input DAG: `align | sort | call` (3-task pipeline)

Output (`function.wdl`):

```wdl
version 1.1

task align {
    input {
        File fastq1
        File fastq2
        String outbase
        File ref
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
    }

    command <<<
        bwa mem -t ~{cpu} ~{ref} ~{fastq1} ~{fastq2} | samtools view -b - > ~{outbase}.bam
    >>>

    output {
        File bam = "~{outbase}.bam"
    }

    requirements {
        cpu: 2
        container: "djhshih/seqkit:0.1"
        time_minutes: 3150
    }
}

workflow main {
    input {
        File fastq1
        File fastq2
        String outbase
        File ref
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_fai
        File ref_pac
        File ref_sa
    }

    call align {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            outbase = outbase,
            ref = ref,
            ...
    }
    call sort {
        input:
            outbase = outbase,
            bam = align.bam,
    }
    call call {
        input:
            outbase = outbase,
            ref = ref,
            bam = sort.bam,
    }

    output {
        File bai = sort.bai
        File bam = sort.bam
        File bcf = call.bcf
    }
}
```

### Mapped step with scatter (`map.swl`)

```wdl
workflow main {
    input {
        String xs
    }

    scatter (call_variant_i in range(length(xs))) {
        call call_variant {
            input:
                fastq1 = fastq1,
                fastq2 = fastq2,
                outbase = outbase,
                ref = ref,
                ...
        }
    }

    output {
        Array[File] bai = call_variant.bai
        Array[File] bam = call_variant.bam
        Array[File] bcf = call_variant.bcf
    }
}
```
