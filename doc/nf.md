# Nextflow DSL2 Transpiler: Implementation Status

## Overview

The Nextflow transpiler (`python/swl/transpile/nf/`) converts compiled SWL DAGs into Nextflow DSL2 `process` + `workflow` blocks. It follows the same architecture as the CWL, WDL, and Snakemake transpilers. The transpiler is ~344 lines of Python (`emit.py`).

**Status:** Implemented. Covers simple processes, pipelines, named ports, resources, output path templates, mapped sub-workflows (scatter), `map_by` (grouped scatter), and sub-workflow inlining. Shell variable interpolation (`${var}`) passes through verbatim to Nextflow's template engine without remapping.

---

## Architecture

```
transpile/nf/
├── __init__.py     # exports transpile_dag_dict, transpile_dag_file
├── __main__.py     # python -m entry point
├── cli.py          # CLI wrapper (uses transpile/_cli.py)
└── emit.py         # ~344 lines, the transpiler
```

### Entry points

- `transpile_dag_file(path)` — read DAG JSON from disk, transpile
- `transpile_dag_dict(data, workflow_id, _top_level)` — main transpilation from parsed DAG data

### Key internal functions

| Function | Purpose |
|----------|---------|
| `_task_to_process(step)` | Convert a task StepCall to a `process { }` block with `input:`, `output:`, `script:`, and directives |
| `_emit_directives(step)` | Map SWL `run` (cpu, memory, time, image) to Nextflow `cpus`, `memory`, `time`, `container` |
| `_interp_to_nf(value)` | Render a SWL `word` interpolation dict into a Nextflow string template `"${var}"` |
| `_dag_to_nf(dag, workflow_id, processes)` | Generate the `workflow { }` block with channel creation, process invocation, and output emission |
| `_binding_to_channel(binding, channels, current_step)` | Convert a DAG binding to a Nextflow channel expression |
| `_input_channel(name, spec)` | Emit `Channel.fromPath()` / `Channel.value()` from workflow input specs |
| `_mapped_step_to_call(step, channels, processes)` | Emit channel wiring for a mapped (scatter) step: join all input channels into a single tuple |
| `_mapped_by_step_to_call(step, channels, processes)` | Emit `groupTuple()` chain for `map_by` grouped scatter |
| `_subworkflow_to_nf(step, parent_id)` | Recursively transpile sub-workflow steps (inline `workflow { }` blocks) |
| `_validate_supported(dag)` | Reject unsupported constructs: literal workflow outputs |

---

## SWL → Nextflow Mapping

### Well Supported

| SWL Concept | Nextflow Output | Details |
|---|---|---|
| **Task step** | `process { }` with `script:` | Script body is embedded literally |
| **Named input ports** | `input:` with `path`/`val` qualifiers | Files → `path name`, strings → `val name`. Mapped steps use `tuple(...)` input |
| **Named output ports** | `output:` with `emit:` name | Output paths use SWL `default` pattern rendered as `"${var}"` template |
| **`run.cpu`** | `cpus <value>` directive | |
| **`run.memory`** | `memory '<value> MB'` directive | Value converted to MB |
| **`run.time`** | `time '<value>m'` directive | Value converted to minutes |
| **`run.image`** | `container '<value>'` directive | |
| **Workflow inputs** | `Channel.fromPath()` / `Channel.value()` | Files → `fromPath(params.name)`, strings → `value(params.name)` |
| **Input binding (workflow input)** | `<name>_ch` channel reference | Channels created from workflow inputs |
| **Input binding (literal)** | `Channel.value(val)` | Inline literal channel |
| **Step output binding** | `PROCESS.out.emit_name` | References named output port |
| **Field projection** | `.map{ it.field }` | Project field from tuple channel |
| **`word` interpolation in output defaults** | Nextflow string template `"${var}"` | `word_interp()` from `common.py` with `var_fn=lambda n: f"${{{n}}}"` |
| **Workflow output** | `<name>_out = <channel>` + `emit:` | Named emit declaration in workflow block |
| **Sub-workflow (no map)** | Inline `workflow { }` block | Recursively transpiled; appears as separate workflow block in output |

### Weakly Supported

| SWL Concept | Status | Issue |
|---|---|---|
| **Mapped step (scatter)** | Single-tuple input | All input channels are `.join()`ed into one tuple, process declares `tuple path(...)` input. Works but consumes all ports in one dimension |
| **`map_by` (grouped scatter)** | `groupKey` + `groupTuple` | Uses Nextflow 22.10+ `groupKey` for keyed grouping. No per-group validation |
| **Mapped sub-workflow outputs** | `.toList()` on output channels | Mapped step outputs use `.toList()` to collect all scatter results into a single channel (may be memory-intensive for many shards) |
| **`[file]` (array of files)** | Falls back to `val` | `to_nf_qualifier()` lacks `[file]` type → emitted as `val` instead of `path` with multiplicity |
| **Sub-workflow (mapped)** | Inline workflow block | Mapped sub-workflows emit a separate `workflow { }` block with duplicated input channels; no module encapsulation |

### Not Supported (explicitly rejected)

| SWL Concept | Behavior |
|---|---|
| **Literal workflow outputs** | Raises `ValueError`: "Nextflow does not support literal workflow outputs" |
| **Record bindings** | Raises `ValueError`: "Record binding with fields [...] must be flattened before Nextflow transpilation" |

### Not Supported (missing implementation)

| SWL Concept | Issue |
|---|---|
| **Shell variable interpolation (`${var}`)** | `${var}` passes through verbatim to Nextflow `script:` block. Nextflow's template engine treats `${var}` as its own expression, so `${cpu}` would fail at runtime (`cpu` is not a Nextflow variable). Should convert to `{input.var}` / `{output.var}` / `${task.cpus}` / `${task.memory}`, similar to the SMK transpiler's `_interpolate_shell()` |
| **`${memory / cpu}` expression** | Passes through verbatim; same interpolation issue as simple `${var}` |
| **Merge binding** | Not rejected in validation but not handled in `_binding_to_channel`. However, Merge bindings are expected to be flattened before DAG generation |

### Nextflow-specific details

- **Channel wiring**: Workflow input channels are created at the top of the `workflow { }` block. Step inputs are wired by joining channels with `.join()`. If a binding is None (unbound), the input name is assumed to be a workflow-level input channel.
- **Mapped processes**: Use a single `tuple(...)` input with all ports joined via `.join()`. This avoids managing multiple scatter dimensions.
- **`map_by`**: Uses `groupKey(value)` + `.groupTuple()` for key-based grouping (Nextflow 22.10+). The grouping is over the source channel elements.
- **Sub-workflows**: Emitted as standalone `workflow { }` blocks in the same output file, with separate channel declarations.
- **Output wiring**: Mapped step outputs are collected via `.toList()` before emission. Non-mapped step outputs are wired directly.
- **Shell interpolation**: The script body is NOT interpolated by the transpiler. `${var}` references pass through literally. In Nextflow DSL2, the `script:` block treats `${var}` as a template expression, so variables in the shell script must match Nextflow process variables (input names, output names, or `${task.*}` for directives).

---

## Support Matrix

| Feature | Status | Tests |
|---------|--------|-------|
| Simple process | ✅ Well supported | function, pipe, explicit |
| Pipeline chain (|) | ✅ Well supported | function, pipe (equivalence) |
| Named I/O ports | ✅ Well supported | All |
| Resources (cpu/memory/time/image) | ✅ Well supported | All |
| Workflow input channels | ✅ Well supported | All |
| Output path templates (`"${var}"`) | ✅ Well supported | function, pipe, explicit |
| Literal input bindings | ✅ Well supported | partial (Literal refs) |
| Field projections | ✅ Well supported | All (channel wiring) |
| Workflow `emit:` declarations | ✅ Well supported | All |
| Sub-workflow (no map) | ✅ Supported | partial, import_partial |
| Mapped step (scatter) | ⚠️ Weak support | map, panel |
| `map_by` (grouped scatter) | ⚠️ Weak support | map_by |
| `[file]` array input type | ❌ Falls to `val` | — |
| Shell variable interpolation (`${var}`) | ❌ Passes through verbatim | — |
| `memory/cpu` expression in shell | ❌ Passes through verbatim | — |
| Literal workflow outputs | ❌ Rejected | — |
| Record bindings | ❌ Rejected | — |

### Coverage

The test suite (`bash test.sh`):
- Generates `tests/nf/*.nf` from all compiled DAGs
- Compares `pipe.nf` and `explicit.nf` against `function.nf` for equivalence of golden output
- `panel.nf`, `map.nf`, `map_by.nf` are generated but have no golden comparison

---

## Known Limitations

1. **Shell variable interpolation**: Nextflow's `script:` block uses `${var}` for its own template expressions. SWL shell scripts with `${cpu}` or `${ref}` pass through verbatim, but Nextflow will try to resolve these as template variables. Since `cpu` is only defined as a `cpus <value>` directive (accessed via `${task.cpus}`), and `ref` is an input (accessed via `${ref}` in the template scope), unqualified `${cpu}` would fail at runtime. A future fix should convert:
   - `${cpu}` → `${task.cpus}`
   - `$var` / `${var}` for inputs → `${var}` (stays same, Nextflow template resolves from input scope)
   - `${memory / cpu}` → `${task.memory / task.cpus}`

2. **`[file]` type qualifier**: Array-of-file types get the default `val` qualifier because `to_nf_qualifier()` in `types.py` only handles `file`, `str`, `int`, `float`. Should emit `path` (or `path` with multiplicity) for `[file]` types.

3. **Mapped sub-workflow output collection**: `.toList()` on mapped step outputs collects ALL results into memory. For large scatter jobs this may cause OOM. A streaming approach (e.g., channel of channels) would be more scalable but is not implemented.

4. **Sub-workflow module isolation**: Sub-workflows are inlined as separate `workflow { }` blocks with duplicated channel declarations. They are not extracted into standalone module files with `include` statements. This works for single-file outputs but doesn't scale to multi-module projects.

5. **Channel join ordering**: `.join()` depends on matching tuple ordering across channels. When inputs come from independent sources (e.g., separate `Channel.fromPath()` globs), element ordering must match. The transpiler does not insert `.sort()` calls to ensure consistent ordering.

---

## Example Output

### Simple pipeline (`function.swl`)

Input DAG: `align | sort | call` (3-task pipeline)

Output (`function.nf`):

```nextflow
process ALIGN {
    cpus 2
    time '3150m'
    container 'djhshih/seqkit:0.1'

    input:
    path fastq1
    path fastq2
    val outbase
    path ref
    path ref_amb
    path ref_ann
    path ref_bwt
    path ref_pac
    path ref_sa

    output:
    path "${outbase}.bam", emit: bam

    script:
    """
    bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam
    """
}

workflow {
    fastq1_ch = Channel.fromPath(params.fastq1, checkIfExists: true)
    fastq2_ch = Channel.fromPath(params.fastq2, checkIfExists: true)
    outbase_ch = Channel.value(params.outbase)
    ref_ch = Channel.fromPath(params.ref, checkIfExists: true)
    ...

    align_ch = fastq1_ch.join(fastq2_ch).join(outbase_ch).join(ref_ch)...
    ALIGN(align_ch)
    sort_ch = ALIGN.out.bam.join(outbase_ch)
    SORT(sort_ch)
    ...

    bam_out = SORT.out.bam
    emit: bam_out
    bcf_out = CALL.out.bcf
    emit: bcf_out
}
```
