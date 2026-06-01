# Nextflow DSL2 Transpiler: Implementation Plan

Transpile SWL compiled DAG JSON to Nextflow DSL2 workflows. Follows the same pattern as `swl/transpile/cwl/` — reads DAG JSON, emits target code, never reads source files.

---

## 1. Nextflow DSL2 Concepts

| Concept | Role | SWL/DAG equivalent |
|---------|------|--------------------|
| `process` | A task with script + directives | `StepCall` with `task.body` |
| `workflow` | Composition block wiring processes | `DAG.inputs` + `DAG.steps` + `DAG.outputs` |
| `take` / `emit` | Module interface (input/output channels) | Workflow inputs/outputs |
| `Channel.of()` | Emit literal values from params | `Input` binding |
| `Channel.fromPath()` | File channel from glob | `file`-typed workflow input |
| `process.out` | Named output channel | `Field(step, output_name)` |
| `.map{ it.field }` | Field projection on tuple/object | `{source: "field", ...}` |
| `.join()` | Merge channels by key | Record merge |
| `.combine()` | Cross product of channels | N/A |
| `.groupTuple()` | Group channel elements by key | `map_by` grouping |
| `.flatten()` | Flatten nested collections | N/A |
| `.toSortedList()` | Collect all elements | Collecting mapped outputs |
| `each` | Parameterize over a list | `map` / scatter |
| `input:` / `output:` | Process port declarations | Step input/output params |
| `val` / `path` / `env` | Input type qualifiers | SWL types (str → val, file → path) |
| `cpus` / `memory` / `time` / `container` | Process directives | `step.run` parameters |
| `publishDir` | Output file publishing | Output path interpolation |
| `workflow.onComplete` | Completion callback | Post-run checks (Sec 7) |
| `include` / `MODULE.nf` | Importing sub-workflows | Workflow step with `definition.dag` |

---

## 2. Transpiler Structure

Follow the CWL module (`swl/transpile/cwl`):

```
swl/transpile/nextflow/
    __init__.py       # package marker
    cli.py            # CLI entry point: nf transpile <dag.json> -o <output.nf>
    __main__.py        # python -m swl.transpile.nextflow
    emit.py            # core transpilation logic
    test_emit.py       # tests
```

### Entry point

```python
# cli.py
import argparse, json, sys
from swl.ir.dag import DAG
from swl.transpile.nextflow.emit import transpile_dag_dict

def main():
    ap = argparse.ArgumentParser('Transpile compiled DAG JSON to Nextflow DSL2')
    ap.add_argument('input', help='compiled DAG json path')
    ap.add_argument('-o', '--output', help='output .nf path')
    args = ap.parse_args()
    data = json.load(open(args.input))
    nf = transpile_dag_dict(data)
    if args.output:
        with open(args.output, 'w') as f:
            f.write(nf)
            f.write('\n')
    else:
        print(nf)
```

### `emit.py` structure

```python
# emit.py — top-level structure
def transpile_dag_dict(data, workflow_id='main'):
    dag = DAG.from_dict(data)
    _validate_supported(dag)
    
    processes = {}
    for step in dag.steps:
        if step.id not in processes:
            processes[step.id] = _step_to_process(step)
    
    workflow = _dag_to_nf(dag, workflow_id, processes)
    
    lines = []
    # 1. Emit process definitions
    for name, body in processes.items():
        lines.append(body)
        lines.append('')
    # 2. Emit workflow block
    lines.append(workflow)
    return '\n'.join(lines)
```

---

## 3. DAG-to-Nextflow Mapping

### 3.1 Task step → Process

A `StepCall` with `type: "task"` becomes a `process` block.

```
process ALIGN {

    // directives from step.run
    cpus 4
    memory '8 GB'
    container 'djhshih/seqkit:0.1'
    publishDir params.outdir, mode: 'copy', pattern: '*.{bam,bai}'

    input:
    path fastq1
    path fastq2
    path ref
    val outbase

    output:
    path "${outbase}.bam", emit: bam

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${fastq1} ${fastq2} |
        samtools view -b - > ${outbase}.bam
    """
}
```

**Mapping rules:**

| DAG field | Nextflow |
|-----------|----------|
| `id` | Process name (upper-cased, sanitized) |
| `run.cpu.value` | `cpus <value>` directive |
| `run.memory.value` | `memory '<value> MB'` directive |
| `run.image.value` | `container '<value>'` directive |
| `run.time.value` | `time '<value>m'` directive |
| `inputs` | `input:` block, each param → `path` or `val` based on type |
| `outputs` | `output:` block, each param → `path` or `val`, with `emit:` name |
| `outputs[].default` | file path in `path` declaration — glob or expression |
| `task.body` | `script:` block (literally embedded) |
| `bindings` | Channel wiring in the `workflow` block (not the process) |

**Type-to-qualifier mapping:**

```python
def _input_qualifier(swl_type):
    return {  # type → (qualifier, nextflow_type)
        'file': ('path', None),
        'str': ('val', 'string'),
        'int': ('val', 'integer'),
        'float': ('val', 'float'),
    }.get(swl_type, ('val', 'string'))
```

**Process name sanitization:**

```python
def _process_name(step_id):
    name = step_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'PROCESS'
    if name[0].isdigit():
        name = '_' + name
    return name.upper()
```

### 3.2 Workflow block

The top-level workflow block wires processes together via channels.

```
workflow {
    // Workflow inputs → channels
    fastq1_ch = Channel.fromPath(params.fastq1, checkIfExists: true)
    ref_ch    = Channel.fromPath(params.ref, checkIfExists: true)
    outbase_ch = Channel.value(params.outbase)

    // Task invocations with implicit channel wiring
    ALIGN(fastq1_ch, fastq2_ch, ref_ch, outbase_ch)
    // ...
}
```

**Channel creation from workflow inputs:**

| Input type | Channel factory |
|------------|----------------|
| `file` (single) | `Channel.fromPath(params.<name>, checkIfExists: true)` |
| `str` | `Channel.value(params.<name>)` |
| `int` | `Channel.value(params.<name>)` |
| `[file]` (array) | `Channel.fromPath(params.<name>, checkIfExists: true).toList()` |

**Binding-to-channel wiring:**

| Binding `source` | Nextflow expression |
|---|---|
| `input` | `<name>_ch` (pre-created workflow input channel) |
| `step`/`step_output` | `<STEP_NAME>.out.<output_name>` |
| `literal` | `Channel.value(<value>)` |
| `field` | `.map{ it.<field_name> }` on the source channel |
| `merge` | `.join()` on the two source channels (matched by ordering) |
| `record` | `.map{ tuple(<fields>) }` combining individual channels |

### 3.3 Channel wiring implementation

The transpiler builds a channel dictionary while emitting the workflow block:

```python
def _dag_to_nf(dag, workflow_id, processes):
    lines = ['workflow {']
    
    # First pass: create channels for workflow inputs
    channels = {}
    for name, spec in dag.inputs.items():
        ch_name = _channel_name(name)
        channels[name] = ch_name
        lines.append(f'    {ch_name} = {_input_channel(name, spec)}')
    
    # Second pass: wire each step's inputs
    for step in dag.steps:
        if getattr(step, 'map', None) is not None:
            lines.append('')
            lines.extend(_mapped_step_to_call(step, channels, processes))
            continue
        
        inputs = []  # ordered list of channel expressions
        for input_name, binding in step.bindings.items():
            ch_expr = _binding_to_channel(binding, channels, step)
            inputs.append(ch_expr)
        
        pname = _process_name(step.id)
        args = ', '.join(inputs)
        lines.append(f'')
        pname = _process_name(step.id)
        
        # Store output channels for downstream consumers
        for out_name in step.outputs:
            channels[f'{step.id}.{out_name}'] = f'{pname}.out.{_nf_emit_name(out_name)}'
    
    # Third pass: wire workflow outputs
    for name, binding in dag.outputs.items():
        ch_expr = _binding_to_channel(binding, channels, None)
        lines.append(f'    {name}_out = {ch_expr}')
        lines.append(f'    emit: {name}_out')
    
    lines.append('}')
    return '\n'.join(lines)
```

**Binding-to-channel expression:**

```python
def _binding_to_channel(binding, channels, current_step):
    """Convert a DAG binding value to a Nextflow channel expression."""
    if isinstance(binding, Input):
        return channels[binding.name]
    
    if isinstance(binding, Literal):
        val = json.dumps(binding.value)
        return f'Channel.value({val})'
    
    if isinstance(binding, Field):
        source_ch = _binding_to_channel(binding.source, channels, current_step)
        if isinstance(binding.source, (StepCall, MappedStep)):
            pname = _process_name(binding.source.id)
            return f'{pname}.out.{_nf_emit_name(binding.name)}'
        return f'{source_ch}.map{{ it.{binding.name} }}'
    
    if isinstance(binding, Merge):
        left = _binding_to_channel(binding.left, channels, current_step)
        right = _binding_to_channel(binding.right, channels, current_step)
        return f'{left}.join({right})'
    
    if isinstance(binding, Record):
        field_exprs = []
        for fname, fbinding in binding.fields.items():
            fch = _binding_to_channel(fbinding, channels, current_step)
            field_exprs.append(f'{fname}: {fch}')
        return f'Channel.of([{", ".join(field_exprs)}])'
    
    raise ValueError(f'Unsupported binding for Nextflow: {type(binding).__name__}')
```

---

## 4. Mapped Steps (`map`) → Scatter

A `MappedStep` without `group_by` represents `map f xs`. In Nextflow this maps to:

### Single-input table (map.source is `input`)

When `map.source` is `{source: "input", name: "xs"}`, the table columns come from individual input arrays. Each column becomes a `Channel.fromPath(...).toList()` or `Channel.value(...).toList()`. The process is wrapped in a `workflow` block that scatters over the zipped columns:

```
workflow {
    fastq1_list = Channel.fromPath(params.fastq1).toList()
    fastq2_list = Channel.fromPath(params.fastq2).toList()
    
    // Zip columns into tuples, scatter with each
    ALIGN(
        fastq1_list.flatMap(),
        fastq2_list.flatMap()
    )
}
```

### Emitted process for mapped step

The mapped process itself is emitted as a regular process (same as a task step), but the calling code uses Nextflow's implicit scatter: when inputs are queues of equal length, each element is processed independently.

For Nextflow DSL2, the simplest scatter is:

```
workflow {
    input_ch = channel.fromPath(params.fastq1)
        .join(channel.fromPath(params.fastq2))
    
    ALIGN(input_ch)
}
```

If `ALIGN` takes `path fastq1, path fastq2` (two inputs), Nextflow must receive them as a single channel of tuples via `.join()`, and the process declares `input: tuple path(fastq1), path(fastq2)`:

```
process ALIGN {
    input:
    tuple path(fastq1), path(fastq2), val(outbase)

    output:
    path "${outbase}.bam", emit: bam

    script:
    """
    bwa mem ... ${fastq1} ${fastq2} > ${outbase}.bam
    """
}
```

**Decision:** When a step has a `map` field, join all its input bindings into a single tuple channel and emit a single input tuple in the process. This avoids having to manage multiple scatter dimensions.

```python
def _mapped_step_to_call(step, channels, processes):
    pname = _process_name(step.id)
    lines = []
    
    # Collect input bindings, join into a single tuple channel
    join_expr = None
    input_names = []
    for input_name in step.inputs:
        binding = step.bindings.get(input_name)
        if binding is None:
            continue  # bound via implicit workflow input
        ch = _binding_to_channel(binding, channels, step)
        input_names.append(input_name)
        if join_expr is None:
            join_expr = ch
        else:
            join_expr = f'{join_expr}.join({ch})'
    
    ch_var = f'{step.id}_ch'
    lines.append(f'    {ch_var} = {join_expr}')
    lines.append(f'    {pname}({ch_var})')
    
    # Store output channels
    for out_name in step.outputs:
        channels[f'{step.id}.{out_name}'] = f'{pname}.out.{_nf_emit_name(out_name)}'
    
    return lines
```

**Process signature for mapped step** — use a single tuple input:

```
process ALIGN {
    input:
    tuple path(fastq1), path(fastq2), val(outbase)
    
    output:
    path "${outbase}.bam", emit: bam
    
    script:
    """
    ...
    """
}
```

---

## 5. `map_by` → `groupTuple` + per-group process

The `map_by` DAG has `map.group_by` set. Nextflow's `groupTuple()` natively supports grouped scatter.

### Strategy: `groupTuple` by key, then process per group

For `map_by align "sample" samples` where `samples` has columns `sample, fastq1, fastq2`:

```
workflow {
    // 1. Create per-row channel of tuples
    rows_ch = Channel
        .fromPath(params.fastq1)
        .join(Channel.fromPath(params.fastq2))
        .join(Channel.value(params.sample))
        .map { fastq1, fastq2, sample ->
            tuple(sample, fastq1, fastq2)
        }
    
    // 2. Group by sample
    groups_ch = rows_ch
        .groupTuple(by: 0)  // group by first element (sample)
        .map { sample, fastq1_list, fastq2_list ->
            tuple(sample, fastq1_list, fastq2_list)
        }
    
    // 3. Process per group
    ALIGN_GROUPED(groups_ch)
}
```

The grouped process `ALIGN_GROUPED` receives arrays (the group) and internally iterates:

```
process ALIGN_GROUPED {
    input:
    tuple val(sample), path(fastq1_list), path(fastq2_list)
    
    output:
    path "*.bam", emit: bam
    
    script:
    """
    # Iterate over group rows
    for i in \$(seq 0 \$((\${#fastq1_list[@]} - 1))); do
        bwa mem ... \${fastq1_list[\$i]} \${fastq2_list[\$i]} > sample_\${sample}_\$i.bam
    done
    """
}
```

### Emitting `map_by`

```python
def _mapped_by_step_to_call(step, channels, processes):
    map_info = step.map
    source = map_info.get('source', {})
    group_key = map_info.get('group_by')
    
    # 1. Build the per-row channel by joining column bindings
    col_names = list(step.input_schema.keys()) if step.input_schema else []
    if group_key not in col_names:
        col_names = [group_key] + col_names
    
    join_expr = None
    for col_name in col_names:
        binding = step.bindings.get(col_name)
        if binding is None and source.get('columns'):
            binding = source['columns'].get(col_name)
        if binding is None:
            continue
        ch = _binding_to_channel(binding, channels, step)
        join_expr = ch if join_expr is None else f'{join_expr}.join({ch})'
    
    ch_var = f'{step.id}_rows'
    lines = [f'    {ch_var} = {join_expr}']
    
    # 2. Group by key column (assume it's the first tuple position)
    key_idx = col_names.index(group_key)
    grouped_var = f'{step.id}_groups'
    lines.append(f'    {grouped_var} = {ch_var}')
    lines.append(f'        .map{{ tuple -> tuple }')  # identity map to ensure tuple
    lines.append(f'        .groupTuple(by: {key_idx})')
    
    # 3. Call process with grouped channel
    pname = _process_name(step.id)
    lines.append(f'    {pname}({grouped_var})')
    
    for out_name in step.outputs:
        channels[f'{step.id}.{out_name}'] = f'{pname}.out.{_nf_emit_name(out_name)}'
    
    return lines
```

---

## 6. Sub-workflow Steps

A step with `type: "workflow"` has a `definition` containing a recursive DAG. In Nextflow DSL2, sub-workflows are emitted as separate modules and `include`d.

```python
def _step_to_process(step):
    if step.type == 'workflow':
        return _subworkflow_to_module(step)
    return _task_to_process(step)
```

**Sub-workflow as a module:**

```
// Generated: _sub_pipeline.nf (or inline)
workflow SUB_PIPLINE {
    take:
    fastq1
    fastq2
    ref
    
    main:
    SORT(fastq1, fastq2, ref)
    CALL(SORT.out.bam, ref)
    
    emit:
    CALL.out.bcf
}
```

**Inclusion in parent workflow:**

```
include { SUB_PIPELINE } from './_sub_pipeline.nf'

workflow {
    SUB_PIPELINE(fastq1_ch, fastq2_ch, ref_ch)
}
```

For simplicity, the initial implementation emits all sub-workflows inline (as nested `workflow` blocks within the same file). Module extraction can be added later.

```python
def _subworkflow_to_nf(step, workflow_id):
    dag = step.task.get('dag', {})
    if not dag:
        dag = step.definition.get('dag', {})
    # recursively transpile
    return transpile_dag_dict(dag, workflow_id=f'{workflow_id}_{step.id}')
```

---

## 7. Resource and Directive Mapping

```python
def _emit_directives(step):
    directives = []
    run = step.run if isinstance(step, StepCall) else getattr(step, 'run', {})
    
    for name, spec in run.items():
        value = spec.get('value') if isinstance(spec, dict) else getattr(spec, 'value', None)
        if value is None:
            continue
        
        if name == 'cpu':
            directives.append(f'    cpus {value}')
        elif name == 'memory':
            directives.append(f"    memory '{value} MB'")
        elif name == 'time':
            directives.append(f"    time '{value}m'")
        elif name == 'image':
            directives.append(f"    container '{value}'")
    
    return '\n'.join(directives)
```

---

## 8. Interpolation Handling

Interpolation words in output defaults and step bindings are converted to Nextflow string templates.

```python
def _interp_to_nf(value):
    """Convert DAG interpolation dict to Nextflow string template."""
    if value is None:
        return None
    if value.get('kind') == 'word':
        parts = value.get('parts', [])
        result = ''
        for part in parts:
            if part.get('kind') == 'literal':
                result += part['text']
            elif part.get('kind') == 'var':
                result += f"${{{part['name']}}}"
            elif part.get('kind') == 'expr':
                result += f"${{{part['text']}}}"
        return result
    return None
```

**Usage in output declarations:**

```python
default_expr = _interp_to_nf(spec.get('default'))
if default_expr:
    output_line = f'path "{default_expr}", emit: {name}'
```

---

## 9. Workflow Output Wiring

Output bindings at the `DAG.outputs` level are converted to `emit:` declarations in the workflow block.

| Output binding | Nextflow emit |
|---|---|
| `{step: "align", output: "bam"}` | `emit: bam` → `ALIGN.out.bam` |
| `{source: "input", name: "fastq1"}` | `emit: fastq1` → params passthrough |
| `{source: "literal", value: 42}` | `emit: result` → `Channel.value(42)` |

```python
def _output_to_emit(name, binding, channels):
    if isinstance(binding, Field) and isinstance(binding.source, (StepCall, MappedStep)):
        pname = _process_name(binding.source.id)
        out_emit = _nf_emit_name(binding.name)
        return f"    {pname}.out.{out_emit} | view{{ it }}"
    return f'    Channel.value(...)'
```

---

## 10. Implementation Order

### Phase 1: Scaffolding + basic task processes

1. Create `swl/transpile/nextflow/` package with `__init__.py`, `cli.py`, `__main__.py`.
2. Implement `transpile_dag_dict()` skeleton: validate, iterate steps, call `_validate_supported`.
3. Implement `_task_to_process()`: emit a Nextflow `process` block from a `StepCall` with all directives, inputs, outputs, and script.
4. Test: round-trip a simple single-task DAG (like `partial.swl` from the CWL tests).

### Phase 2: Workflow block + channel wiring

5. Implement `_dag_to_nf()`: emit the top-level `workflow { }` block.
6. Implement `_binding_to_channel()`: handle `input`, `literal`, `step_output`, `field` bindings.
7. Implement `_input_channel()`: emit `Channel.fromPath()` / `Channel.value()` from `InputSpec`.
8. Implement workflow output wiring.
9. Test: transpile `function.swl` (three tasks in a pipeline) and verify channel wiring.

### Phase 3: Mapped steps (map)

10. Implement `_mapped_step_to_call()`: join column channels into tuples, emit process call.
11. Modify `_task_to_process()` to emit `tuple` input when step has a `map` field.
12. Test: transpile `batch.swl` (map align over samples) and verify scatter.

### Phase 4: map_by (grouped scatter)

13. Implement `_mapped_by_step_to_call()`: emit `.groupTuple(by: ...)` chain.
14. Test: transpile a simple `map_by` workflow and verify `groupTuple` call.

### Phase 5: Sub-workflows

15. Implement `_subworkflow_to_module()`: recursively transpile `definition.dag`.
16. Implement `include` statement emission for sub-workflow references.
17. Test: transpile `import_partial.swl` (workflow step containing sub-workflow).

### Phase 6: Interpolation + run parameters

18. Implement `_interp_to_nf()` for output path templates.
19. Wire interpolation into process `output:` declarations.
20. Verify `test_output_glob_uses_cwl_expression` equivalent works for Nextflow.

### Phase 7: Record merge + field projection

21. Implement `Merge` binding → `.join()` emission.
22. Implement `Record` binding → tuple construction.
23. Implement nested field projection → `.map{ it.field }`.

### Phase 8: Edge cases and hardening

24. Handle special characters in process names.
25. Handle empty inputs.
26. Handle `map_by` with lambda (generated sub-workflow).
27. Handle partial application root `map_by`.
28. Add comprehensive validation in `_validate_supported()`.

---

## 11. Validation (`_validate_supported`)

```python
def _validate_supported(dag):
    for step in dag.steps:
        if getattr(step, 'map', None) is not None:
            map_info = step.map
            source = map_info.get('source', {})
            group_by = map_info.get('group_by')
            
            if group_by is not None:
                # group_by is supported natively via groupTuple
                if source.get('source') != 'input' and source.get('source') != 'table':
                    raise ValueError(f'Nextflow map_by requires input or table source, got {source.get("source")}: {step.id}')
                if group_by not in (step.input_schema or {}):
                    raise ValueError(f'Group key {group_by} not in input schema for {step.id}')
        
        for name, binding in step.bindings.items():
            _validate_binding(binding, step)
        
        for name, binding in getattr(dag, 'outputs', {}).items():
            _validate_output_binding(binding)


def _validate_binding(value, step):
    if isinstance(value, Merge):
        raise ValueError(f'Nextflow does not support merge bindings (step {step.id})')
    if isinstance(value, Record):
        # Records in bindings are OK if they are simple field-value pairs
        pass
```

---

## 12. Feature Support Summary

| DAG feature | Nextflow support | Status |
|---|---|---|
| Task step | `process` with `script:` | Planned Phase 1 |
| Input binding | `Channel.fromPath()` / `Channel.value()` | Planned Phase 2 |
| Step output binding | `PROCESS.out.emit_name` | Planned Phase 2 |
| Literal binding | `Channel.value(val)` | Planned Phase 2 |
| Field projection | `.map{ it.field }` | Planned Phase 7 |
| Record merge | `.join()` | Planned Phase 7 |
| Record binding | Tuple construction | Planned Phase 7 |
| `map` (scatter) | `.join()` → tuple input with implicit scatter | Planned Phase 3 |
| `map_by` | `.groupTuple(by: idx)` | Planned Phase 4 |
| Sub-workflow | `include` + `workflow` block | Planned Phase 5 |
| CPU directive | `cpus <value>` | Planned Phase 1 |
| Memory directive | `memory '<val> MB'` | Planned Phase 1 |
| Docker image | `container '<val>'` | Planned Phase 1 |
| Time directive | `time '<val>m'` | Planned Phase 1 |
| Output glob | `path "${var}.bam"` | Planned Phase 6 |
| Interpolation (var) | `"${var}"` in path templates | Planned Phase 6 |
| Interpolation (expr) | `"${expr}"` pass-through | Planned Phase 6 |
| Table binding | Channel-of-tuples from column joins | Planned Phase 3 |

---

## 13. Key Risks and Mitigations

| Risk | Mitigation |
|---|---|
| **Channel ordering**: Nextflow channels are ordered queues, but `.join()` requires matching order. If input channels have different orderings (e.g., from file glob), elements may mismatch. | Use `.join()` only on channels derived from the same source. For independent file globs, sort both by filename: `Channel.fromPath(params.fastq1).sort()`. |
| **File staging**: Nextflow stages input files into the process work directory. The DAG JSON stores abstract paths; the transpiler must not resolve them. | Use `path` qualifier (not `val` with path string). Nextflow handles staging automatically. |
| **`groupTuple` memory**: All rows for a group are buffered in memory. Large groups may cause OOM. | Document this limitation. Recommend pre-grouping by the workflow author for large datasets. |
| **Process name collisions**: Multiple steps may generate the same sanitized name. | Append numeric suffix for duplicates: `ALIGN`, `ALIGN_1`, `ALIGN_2`. |
| **Script interpolation conflicts**: Nextflow uses `$` and `${}` for its own string templates. SWL bash scripts also use `${}`. | Escape `$` in process `script:` blocks when it represents shell variable, not Nextflow variable. Use `script:` (not `shell:`) and let Nextflow pass through unescaped `$`. |
| **Sub-workflow parameter scope**: Sub-workflows need their own `take`/`emit` interface. Parameters from the outer workflow don't leak in. | Always generate explicit `take:` and `emit:` declarations for sub-workflow modules. |
| **`map_by` with grouped tuple size**: `groupTuple` by index assumes the key is at a fixed position in the tuple. If the schema changes (e.g., optional input), the index breaks. | Use named `.groupTuple(by: <name>)` in DSL2 (requires Nextflow 22.10+). Fall back to index-based for older versions. |

---

## 14. Example Output: Complete Workflow

### Input DAG (`function.swl` compiled)

Simple pipeline: `align | sort | call`

### Output Nextflow (`function.nf`)

```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastq1 = ''
params.fastq2 = ''
params.ref = ''
params.ref_fai = ''
params.outbase = 'output'

process ALIGN {
    cpus 2
    memory '8192 MB'
    container 'djhshih/seqkit:0.1'
    time '3000m'

    input:
    path fastq1
    path fastq2
    path ref
    path ref_fai
    val outbase

    output:
    path "${outbase}.bam", emit: bam
    path "${outbase}.bai", emit: bai

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${fastq1} ${fastq2} |
        samtools view -b - > ${outbase}.bam
    samtools index ${outbase}.bam ${outbase}.bai
    """
}

process SORT {
    cpus 2
    memory '8192 MB'
    container 'djhshih/seqkit:0.1'

    input:
    path bam
    val outbase

    output:
    path "${outbase}.sorted.bam", emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -m 8G ${bam} -o ${outbase}.sorted.bam
    """
}

process CALL {
    container 'djhshih/seqkit:0.1'

    input:
    path bam
    path ref

    output:
    path "*.bcf", emit: bcf

    script:
    """
    bcftools mpileup -Ou -f ${ref} ${bam} |
        bcftools call -mv -Ob -o variants.bcf
    """
}


workflow {
    // Workflow inputs
    fastq1_ch = Channel.fromPath(params.fastq1, checkIfExists: true)
    fastq2_ch = Channel.fromPath(params.fastq2, checkIfExists: true)
    ref_ch    = Channel.fromPath(params.ref, checkIfExists: true)
    ref_fai_ch = Channel.fromPath(params.ref_fai, checkIfExists: true).optional()
    outbase_ch = Channel.value(params.outbase)

    // Pipeline: align | sort | call
    // align
    align_input = fastq1_ch.join(fastq2_ch).join(ref_ch).join(ref_fai_ch).join(outbase_ch)
    ALIGN(align_input)

    // sort (x // a): merge workflow input with align output
    sort_input = fastq1_ch.join(fastq2_ch).join(ref_ch).join(ref_fai_ch).join(outbase_ch)
        .join(ALIGN.out.bam)
        .map { fastq1, fastq2, ref, ref_fai, outbase, bam ->
            tuple(outbase, bam)
        }
    SORT(sort_input)
}
```
