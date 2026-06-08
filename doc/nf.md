# Nextflow DSL2 Transpiler: Implementation Status

## Overview

The Nextflow transpiler (`python/swl/transpile/nf/`) converts compiled SWL DAGs into Nextflow DSL2 `process` + `workflow` blocks. The transpiler is ~519 lines of Python (`emit.py`).

**Status:** Implemented. Covers simple processes, pipelines, named ports, resources, output path templates, shell variable interpolation with scope resolution (inputs → `${var}`, run params → `${task.cpus}`/`${task.memory}`/`${task.time}`), sub-workflow inlining (mapped and unmapped), per-field channel wiring, `map_by` (grouped scatter), and input-sourced maps.

---

## Architecture

```
transpile/nf/
├── __init__.py     # exports transpile_dag_dict, transpile_dag_file
├── __main__.py     # python -m entry point
├── cli.py          # CLI wrapper (uses transpile/_cli.py)
└── emit.py         # transpiler
```

### Entry points

- `transpile_dag_file(path)` — read DAG JSON from disk, transpile
- `transpile_dag_dict(data, workflow_id, _top_level)` — main transpilation from parsed DAG data

### Design

- **Variable reference forms**: Both `$var` (unbraced) and `${var}` (braced) are handled identically via `interp_script()`'s regex — they are equivalent bash syntax for simple variable names. `${expr}` (non-word content inside braces) is treated as an expression.
- **Scope resolution**: Input vars keep `${var}` form (Nextflow's `script:` template engine resolves from `input:` scope). Run vars are remapped: `cpu` → `${task.cpus}`, `memory` → `${task.memory}`, `time` → `${task.time}` via `_NF_RUN_MAP`. Shell builtins and unknown vars pass through verbatim.
- **Shell built-in allowlist**: `_BUILTIN_VARS` in `bash.py` prevents shell builtins from being incorrectly resolved as SWL variables.

### Key internal functions

| Function | Purpose |
|----------|---------|
| `_task_to_process(step)` | Convert a task StepCall to a `process { }` block with `input:`, `output:`, `script:`, and directives |
| `_emit_directives(step)` | Map SWL `run` (cpu, memory, time, image) to Nextflow `cpus`, `memory`, `time`, `container` |
| `_interp_to_nf(value)` | Render a SWL `word` interpolation dict into a Nextflow string template `"${var}"` |
| `_collect_processes(dag, processes)` | Recursively collect process names from all nested DAGs |
| `_emit_processes(dag, lines, processes)` | Recursively emit process blocks from all nested DAGs |
| `_dag_to_nf(dag, workflow_id, processes)` | Generate the `workflow { }` block with channel creation, process invocation, and output emission |
| `_inline_dag_steps(dag, channels, lines, processes)` | Process all steps in a DAG, emitting process calls directly (handles mapped, unmapped, and nested workflow steps) |
| `_inline_workflow_step(step, channels, lines, processes)` | Inline a non-mapped workflow step's inner DAG into the current workflow |
| `_binding_to_channel(binding, channels, current_step)` | Convert a DAG binding to a Nextflow channel expression |
| `_input_channel(name, spec)` | Emit `Channel.fromPath()` / `Channel.value()` from workflow input specs |
| `_mapped_step_to_call(step, channels, processes)` | Emit channel wiring for a mapped (scatter) task step |
| `_mapped_workflow_step_to_call(step, channels, processes)` | Emit per-field channels and inline inner processes for a mapped workflow step |
| `_mapped_by_step_to_call(step, channels, processes)` | Emit `groupTuple()` chain for `map_by` grouped scatter |
| `_validate_supported(dag)` | Reject unsupported constructs: literal workflow outputs |

---

## SWL → Nextflow Mapping

### Well Supported

| SWL Concept | Nextflow Output | Details |
|---|---|---|
| **Task step** | `process { }` with `script:` | Script body interpolated via `interp_script()` |
| **Shell variable interpolation (`${var}`)** | `${var}` / `${task.*}` in `script:` | Input vars → `${var}` (Nextflow resolves from `input:`), run params → `${task.cpus}` / `${task.memory}` / `${task.time}` |
| **Shell expression interpolation (`${expr}`)** | `${...}` with resolved run vars | `${memory / cpu}` → `${task.memory / task.cpus}`. Input vars pass through verbatim |
| **Named input ports** | `input:` with `path`/`val` qualifiers | Files → `path name`, strings → `val name`. Mapped tasks use `tuple(...)` input |
| **Named output ports** | `output:` with `emit:` name | Output paths use SWL `default` pattern rendered as `"${var}"` template |
| **`run.cpu`** | `cpus <value>` directive | |
| **`run.memory`** | `memory '<value> MB'` directive | Value converted to MB |
| **`run.time`** | `time '<value>m'` directive | Value converted to minutes |
| **`run.image`** | `container '<value>'` directive | |
| **Workflow inputs** | `Channel.fromPath()` / `Channel.value()` | Files → `fromPath(params.name)`, strings → `value(params.name)` |
| **Input-sourced map inputs** | `Channel.fromList(params.name)` | Array inputs for mapped steps use `fromList` |
| **Input binding (workflow input)** | `<name>_ch` channel reference | Channels created from workflow inputs |
| **Input binding (literal)** | `Channel.value(val)` | Inline literal channel |
| **Step output binding** | `PROCESS.out.emit_name` | References named output port |
| **Field projection** | `.map{ it.field }` | Project field from tuple channel |
| **`word` interpolation in output defaults** | Nextflow string template `"${var}"` | `word_interp()` from `common.py` with `var_fn=lambda n: f"${{{n}}}"` |
| **Shell variable interpolation** | `interp_script()` via `_interpolate_shell` | Input vars → `${var}`, run vars → `${task.cpus}`/`${task.memory}`/`${task.time}`. Expressions resolve inner run vars |
| **Workflow output** | `<name> = <channel>` + `emit:` | Named emit declaration in workflow block (no `_out` suffix) |
| **Sub-workflow (no map)** | Inlined in entry workflow | Inner processes emitted at top level, wired directly in entry workflow; no separate `workflow { }` block |
| **Mapped workflow step** | Per-field channels + inlined pipeline | Struct fields extracted via `.map{ x -> x.field }`, inner processes wired directly as a pipeline; no sub-workflow wrapper |
| **Mapped task step (table-source)** | `.join()` into tuple + process call | Input channels from table columns joined into a single tuple channel for the mapped process |
| **`[file]` type qualifier** | `path` qualifier | `to_nf_qualifier()` now handles `[file]` types |

### Weakly Supported

| SWL Concept | Status | Issue |
|---|---|---|
| **`map_by` (grouped scatter)** | `groupKey` + `groupTuple` | Uses Nextflow 22.10+ `groupKey` for keyed grouping. No per-group validation |
| **Mapped step output collection** | `.toList()` on output channels | Mapped step outputs use `.toList()` to collect all scatter results into a single channel (may be memory-intensive for many shards) |

### Not Supported (explicitly rejected)

| SWL Concept | Behavior |
|---|---|
| **Literal workflow outputs** | Raises `ValueError`: "Nextflow does not support literal workflow outputs" |
| **Record bindings** | Raises `ValueError`: "Record binding with fields [...] must be flattened before Nextflow transpilation" |

### Not Supported (missing implementation)

| SWL Concept | Issue |
|---|---|
| **Merge binding** | Not rejected in validation but not handled in `_binding_to_channel`. However, Merge bindings are expected to be flattened before DAG generation |

### Nextflow-specific details

- **Channel wiring**: Workflow input channels are created at the top of the `workflow { }` block. Input-sourced map sources use `Channel.fromList(params.x)`. Process calls receive multiple channels directly (no `.join()`) — NF pairs them automatically. `.join()` is only used for combining process output channels with data channels (e.g., `SORT(ALIGN.out.bam, outbase)`).
- **Sub-workflows**: Both mapped and unmapped sub-workflows are inlined. Their inner processes are emitted at the top level alongside the entry workflow's direct processes. The entry workflow wires all processes directly — no separate `workflow { }` wrapper blocks are emitted.
- **Mapped workflow steps**: Struct fields are extracted into per-field channels via `.map{ x -> x.field }`. Inner pipeline processes receive these channels directly. Mapped workflow outputs use `.toList()` to gather scatter results.
- **`map_by`**: Uses `groupKey(value)` + `.groupTuple()` for key-based grouping (Nextflow 22.10+). The grouping is over the source channel elements.
- **Output wiring**: Mapped step outputs are collected via `.toList()` before emission. Non-mapped step outputs are wired directly.
- **Shell interpolation**: The script body is interpolated by `_interpolate_shell()` via `interp_script()` from `common.py`. Both `$var` and `${var}` are recognized. Input variables stay as `${var}` (Nextflow resolves from `input:`). Run params (`cpu`, `memory`, `time`) are remapped to `${task.cpus}`, `${task.memory}`, `${task.time}`. Inner variables in `${expr}` expressions are resolved per scope. Shell builtins (`$HOME`) pass through verbatim.

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
| Sub-workflow (no map) | ✅ Inlined | partial, import_partial |
| Mapped workflow step | ✅ Inlined pipeline | map_root, batch_workflow |
| Mapped task step (scatter) | ✅ Supported | batch, batch_lambda |
| `map_by` (grouped scatter) | ⚠️ Weak support | map_by |
| `[file]` array input type | ✅ Supported | — |
| Shell variable interpolation (`${var}`, `$var`) | ✅ Well supported | function (bwa mem line) |
| `memory/cpu` expression in shell | ✅ Supported | function (`${task.memory / task.cpus}`) |
| Literal workflow outputs | ❌ Rejected | — |
| Record bindings | ❌ Rejected | — |

### Coverage

The test suite (`bash test.sh`):
- Generates `tests/nf/*.nf` from all compiled DAGs
- Compares `pipe.nf` and `explicit.nf` against `function.nf` for equivalence of golden output
- `panel.nf`, `map.nf`, `map_by.nf` are generated but have no golden comparison

---

## Known Limitations

1. **Mapped step output collection**: `.toList()` on mapped step outputs collects ALL results into memory. For large scatter jobs this may cause OOM. A streaming approach (e.g., channel of channels) would be more scalable but is not implemented.

2. **Channel join ordering**: `.join()` depends on matching tuple ordering across channels. When inputs come from independent sources (e.g., separate `Channel.fromPath()` globs), element ordering must match. The transpiler does not insert `.sort()` calls to ensure consistent ordering.

3. **Module isolation**: All processes and wiring are emitted into a single file. There is no module extraction with `include` statements. This works for single-file outputs but doesn't scale to multi-module projects.

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
    path ref
    path ref_amb
    path ref_ann
    path ref_bwt
    path ref_pac
    path ref_sa
    val outbase

    output:
    path "${outbase}.bam", emit: bam

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam
    """
}

workflow {
    fastq1_ch = Channel.fromPath(params.fastq1, checkIfExists: true)
    fastq2_ch = Channel.fromPath(params.fastq2, checkIfExists: true)
    outbase_ch = Channel.value(params.outbase)
    ref_ch = Channel.fromPath(params.ref, checkIfExists: true)
    ...

    ALIGN(fastq1_ch, fastq2_ch, ref_ch, ...)
    SORT(ALIGN.out.bam, outbase_ch)
    CALL(SORT.out.bam, ref_ch, ref_fai_ch, outbase_ch)

    bam = SORT.out.bam
    emit: bam
    bcf = CALL.out.bcf
    emit: bcf
}
```

### Mapped pipeline (`map_root.swl`)

Input DAG: `map call_variant` where `call_variant = align | sort | call`

Output:

```nextflow
process ALIGN { ... }
process SORT { ... }
process CALL { ... }

workflow {
    xs_ch = Channel.fromList(params.xs)

    fastq1 = xs_ch.map{ x -> file(x.fastq1) }
    fastq2 = xs_ch.map{ x -> file(x.fastq2) }
    outbase = xs_ch.map{ x -> x.outbase }
    ref = xs_ch.map{ x -> file(x.ref) }
    ref_fai = xs_ch.map{ x -> file(x.ref_fai) }
    ...

    ALIGN(fastq1, fastq2, ref, ...)
    SORT(ALIGN.out.bam, outbase)
    CALL(SORT.out.bam, ref, ref_fai, outbase)

    bam = SORT.out.bam.toList()
    emit: bam
    bai = SORT.out.bai.toList()
    emit: bai
    bcf = CALL.out.bcf.toList()
    emit: bcf
}
```
