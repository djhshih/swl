# Feasibility Report: Snakemake Transpiler for SWL

## Summary

Implementing a Snakemake transpiler for SWL is **highly feasible** with moderate effort (~400–600 lines, comparable to the existing Nextflow transpiler). The main architectural challenge is bridging SWL's **port-based** dependency model (explicit named inputs/outputs connecting steps) with Snakemake's **filename-based** dependency model (implicit resolution via filesystem wildcards).

---

## Snakemake's Execution Model

Snakemake workflows are defined in a `Snakefile` (valid Python extended with rule syntax):

```
rule <name>:
    input:
        "path/to/{wildcard}.txt",
        named_input="other/{wildcard}.bam"
    output:
        "path/to/{wildcard}.out.txt"
    params:
        threshold=0.5
    threads: 8
    resources:
        mem_mb=8000
    container: "docker://ubuntu:latest"
    shell:
        "some_tool {input.named_input} --threshold {params.threshold} > {output}"
```

Key characteristics:
- **Wildcards** (`{sample}`, `{dataset}`) are inferred from output filenames and propagated to inputs automatically — Snakemake builds the DAG by matching output files to input file patterns.
- **Dependencies are implicit**: if rule B's input pattern matches rule A's output pattern, Snakemake links them.
- **Named inputs/outputs** allow port-like access: `{input.bam}`, `{output.vcf}`.
- **Resources** (`threads`, `resources`, `container`) map directly to SWL `run` parameters.
- **Scatter** is implicit via wildcards — no explicit scatter construct needed.
- **`params`** section holds non-file parameters.
- **`expand()`** helper generates multiple filenames from wildcard combinations for aggregation steps.

---

## Mapping SWL Concepts to Snakemake

### 1. Simple Task Rules (high confidence, ~60 lines)

Each SWL step with a `StepCall` maps to a Snakemake `rule` with `shell:` directive:

**SWL DAG step:**
```json
{
  "id": "align",
  "path": "tasks/align.sh",
  "bindings": {"fastq1": {"source": "input", "name": "sample_fastq"}, "ref": {"source": "input", "name": "genome"}},
  "outputs": ["bam"],
  "task": {
    "inputs": {"fastq1": {"type": "file"}, "ref": {"type": "file"}},
    "outputs": {"bam": {"type": "file", "default": {"kind": "word", "parts": [{"kind": "literal", "text": "{sample}.bam"}]}}},
    "run": {"cpu": {"value": 8}, "memory": {"value": 8000}},
    "body": "bwa mem $ref $fastq1 > $bam"
  }
}
```

**Snakemake rule:**
```python
rule align:
    input:
        fastq1=lambda wildcards: f"data/{wildcards.sample}.fastq",
        ref="data/genome.fa"
    output:
        bam="{sample}.bam"
    threads: 8
    resources:
        mem_mb=8000
    shell:
        "bwa mem {input.ref} {input.fastq1} > {output.bam}"
```

The challenge is that SWL uses explicit port wiring, while Snakemake expects filename patterns. Where SWL steps pass entire records, Snakemake needs explicit file paths. This requires the transpiler to generate filename patterns that align across rules — essentially synthesizing wildcard names from SWL's map metadata.

### 2. Pipeline / Chain Desugaring (high confidence, ~30 lines)

SWL pipelines like `align | sort` are already desugared by the compiler into explicit record-passing steps. The transpiler sees explicit step bindings. Each step becomes a rule; Snakemake's DAG resolution links them if output filenames of one match input filenames of the next.

The transpiler needs to choose output paths consistently. Strategy: use `step_id/{wildcard}.{output_name}` patterns so downstream steps can reference `step_id/{wildcard}.{output_name}` as their inputs.

### 3. Named Input/Output Ports (high confidence, ~40 lines)

SWL steps reference inputs by name. Snakemake supports named inputs and outputs natively:

```
rule align:
    input:
        fastq1="data/{sample}.fastq",
        ref="data/genome.fa"
    output:
        bam="results/{sample}.bam"
    shell:
        "bwa mem {input.ref} {input.fastq1} > {output.bam}"
```

Named outputs are accessed via `{output.bam}`. This maps cleanly.

### 4. Resource Requirements (high confidence, ~20 lines)

| SWL run param | Snakemake equivalent |
|---------------|----------------------|
| `cpu` | `threads: <value>` |
| `memory` | `resources: mem_mb=<value>` |
| `time` | `resources: runtime_minutes=<value>` |
| `image` | `container: "docker://<value>"` |

Direct mapping. No structural challenges.

### 5. Map / Scatter (high confidence, ~50 lines)

SWL's `map f xs` scatters function `f` over table `xs`. Snakemake handles this natively through wildcards.

**SWL:** `map align reads` where `reads` is a table with columns `{sample, fastq1, fastq2}`.

**Snakemake:** The `reads` table becomes input file patterns with wildcards:

```python
rule align:
    input:
        fastq1="data/{sample}.R1.fastq",
        fastq2="data/{sample}.R2.fastq"
    output:
        bam="results/{sample}.bam"
    shell:
        "bwa mem {input.fastq1} {input.fastq2} > {output.bam}"
```

For SWL `map` where the table source is an input (not filesystem files), use Snakemake input functions or `params` to pass data. For file-backed tables, wildcards work directly.

### 6. Map By / Grouped Scatter (medium confidence, ~80 lines)

SWL's `map_by f key xs` groups rows by `key` and applies `f` to each group.

Snakemake handles this in two phases:
1. **Grouping rule**: Aggregates per key using `expand()` or `collect` into per-key files.
2. **Processing rule**: Uses the grouping key as a wildcard.

```python
rule group_by_key:
    input:
        data=expand("data/{sample}.fastq", sample=SAMPLES)
    output:
        groups="{key}.txt"
    run:
        # Python code to group by key and write groups
        ...

rule process_group:
    input:
        group="groups/{key}.txt"
    output:
        result="results/{key}.out"
    shell:
        "tool {input.group} > {output}"
```

For `map_by` with table inputs, use Snakemake's `group_by()` function (available since Snakemake 7.x) or implement grouping via Python in a `run:` block. This is the most complex mapping — similar complexity to the CWL and WDL map_by implementations (~80 lines).

### 7. Workflow Inputs (high confidence, ~20 lines)

SWL workflow inputs map to Snakemake config or `params`:

```python
# SWL inputs become config or params
configfile: "config.yaml"

rule all:
    input:
        "results/final.vcf"
```

File-typed inputs can be referenced directly. Non-file inputs (strings, ints) use `params:`.

### 8. Record Bindings (low-medium confidence, ~50 lines)

SWL record literals used in step bindings have no direct Snakemake equivalent. Strategy: flatten record bindings into individual named inputs at the rule level, similar to what the CWL transpiler does with `ExpressionTool`. For Snakemake, this means:

```python
rule example:
    input:
        field1="path/to/field1.txt",
        field2="path/to/field2.txt"
    output:
        result="output.txt"
    shell:
        "tool {input.field1} {input.field2} > {output}"
```

If a record is a workflow output, use a dummy aggregation rule.

### 9. Sub-workflows (medium confidence, ~40 lines)

SWL sub-workflow imports become Snakemake `module` directives:

```python
module other_workflow:
    snakefile:
        "path/to/subworkflow.smk"
    config:
        workflow: "sub"
```

Or inline the sub-workflow's rules with prefixed rule names. The existing WDL and CWL transpilers show the pattern for recursive sub-workflow materialization.

### 10. String Interpolation (high confidence, ~20 lines)

SWL `${var}` interpolation maps to Snakemake `{...}` syntax (Python format minilanguage):

| SWL | Snakemake |
|-----|-----------|
| `$var` or `${var}` | `{params.var}` or `{input.var}` |
| `${expr}` | `{params.expr}` (via `params` block) |

The existing `word_interp()` utility in `transpile/common.py` can be reused with a Snakemake-specific rendering function.

---

## Implementation Plan

### File structure (following existing pattern)

```
transpile/smk/
├── __init__.py          # exports transpile_dag_dict, transpile_dag_file
├── __main__.py          # python -m entry point
├── cli.py               # CLI main, boilerplate
└── emit.py              # ~500 lines, the transpiler
```

### Reusable utilities from `transpile/common.py`

All transpilers share `swl.transpile.common`. The Snakemake transpiler should use these:

| Utility | Usage in Smk transpiler |
|---------|-------------------------|
| `step_name(id, 'rule')` | Sanitize step IDs into rule names: `align` → `align`, `Align_Bwa` → `align_bwa` |
| `workflow_name(id, 'workflow')` | Sanitize workflow name for the `module` directive |
| `emit_name(name)` | Replace `-` with `_` in output/input port names |
| `word_interp(value, lit_fn, var_fn, expr_fn)` | Render SWL interpolation words into `{...}` syntax for `shell:` and output glob strings |
| `run_value(run, 'cpu'\|'memory'\|'time'\|'image')` | Extract resource values for `threads:`, `resources:`, `container:` |
| `source_kind(source)` | Distinguish table vs input sources for `map`/`map_by` |
| `source_input_name(source)` | Get the input name from a map source (for wildcard generation) |
| `table_columns(source)` | Get column names from a table source |
| `column_input_name(source, col)` | Get the input name for a table column |
| `field_chain_parts(value)` | Resolve `Input → Field → Field` chains into root + path |
| `field_path_after_first(value)` | Get the tail segment after the first field access |

### Module: `__init__.py` (1 line)

```python
from .emit import transpile_dag_dict, transpile_dag_file
```

### Module: `__main__.py` (5 lines)

```python
from swl.transpile.smk.cli import main

if __name__ == '__main__':
    main()
```

### Module: `cli.py` (6 lines)

Fully generic, uses `transpile/_cli.py`:

```python
from swl.transpile._cli import run
from swl.transpile.smk.emit import transpile_dag_file


def main():
    run('Transpile compiled DAG JSON to Snakemake', 'smk', transpile_dag_file, False)
```

### Module: `emit.py` — Detailed function breakdown

#### `transpile_dag_file(path) → str` (~3 lines)

Read DAG JSON from disk and delegate:

```python
def transpile_dag_file(path):
    dag = DAG.read(path)
    return transpile_dag_dict(dag.to_dict())
```

#### `transpile_dag_dict(data, workflow_id='main', _top_level=True) → str` (~30 lines)

Main entry: deserialize, validate, collect rules, assemble Snakefile.

```python
from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.transpile.common import (
    column_input_name, emit_name, field_chain_parts, field_path_after_first,
    run_value, source_input_name, source_kind, step_name,
    table_columns, word_interp, workflow_name,
)
from swl.types import to_array_type
```

**Algorithm:**
1. `DAG.from_dict(data)`, call `dag.validate()`.
2. Reject literal workflow outputs (same constraint as WDL & CWL).
3. First pass: collect `StepCall` → rule name mapping. Recursively process sub-workflows (steps with `type == 'workflow'`).
4. Build lines:
   - Rule definitions from task steps (`type != 'workflow'`), each via `_task_to_rule()`.
   - Sub-workflow Snakefiles via `_subworkflow_to_smk()`.
   - Workflow body via `_dag_to_smk()` (an `rule all:` target + wiring).
5. Return `'\n'.join(lines)`.

```
def transpile_dag_dict(data, workflow_id='main', _top_level=True):
    dag = DAG.from_dict(data)
    dag.validate()
    for name, output in dag.outputs.items():
        if isinstance(output.value if isinstance(output, OutputSpec) else output, Literal):
            raise ValueError(f'Snakemake does not support literal workflow outputs: {name}')

    rule_names = {}
    sub_workflows = []
    for step in dag.steps:
        if step.id not in rule_names:
            if step.type == 'workflow':
                smk = _subworkflow_to_smk(step, workflow_id)
                rule_names[step.id] = _rule_name(step.id)
                sub_workflows.append(smk)
            else:
                rule_names[step.id] = _rule_name(step.id)

    lines = []
    for step in dag.steps:
        if step.type == 'workflow':
            continue
        if all(s.type == 'workflow' for s in dag.steps if s.id != step.id):
            continue
        lines.append(_task_to_rule(step, dag, rule_names))
        lines.append('')
    for sw in sub_workflows:
        lines.append(sw)
        lines.append('')
    lines.append(_dag_to_smk(dag, workflow_id, rule_names))
    return '\n'.join(lines)
```

#### `_rule_name(step_id) → str` (~3 lines)

```python
def _rule_name(step_id):
    return step_name(step_id, 'rule')
```

#### `_wf_name(workflow_id) → str` (~3 lines)

```python
def _wf_name(workflow_id):
    return workflow_name(workflow_id, 'workflow')
```

#### `_task_to_rule(step, dag, rule_names) → str` (~80 lines)

Convert a task `StepCall` into a `rule {name}:` block.

**Input handling:**
- Iterate over `task['inputs']` keys.
- For each input name, look up `step.bindings[name]`.
- Dispatch on binding type:
  - `Input` → `input_name="path/to/{wildcard}.ext"` (wildcard if mapped, else literal path).
  - `Field(root=Input)` → chain through field path to build path expression.
  - `Field(root=StepCall)` → reference previous step's output file.
  - `Literal` → if value is a file path, use it; otherwise emit as `params`.
  - `Record` → flatten into individual named inputs.
  - Default → emit as `params` section entry.

**Output handling:**
- Iterate over `step.outputs`.
- Build output paths using the step's `task['outputs'][name]['default']` interpolation via `_interp_to_smk()`.
- If no default, generate `"results/{step_id}/{name}"`.

**Shell command:**
- Use `task.get('body', '')` directly.
- Interpolate bash variables using `_interp_to_smk()` via `params` for non-file values.
- Replace references like `$input_name` with `{input.input_name}` and `$output_name` with `{output.output_name}`.

**Resources:**
- Delegate to `_emit_resources()`.

**Scatter detection:**
- If `step.map` is present and no `group_by`, this step uses wildcards. Add wildcards to output filenames and propagate to inputs.
- If `step.map.group_by` is present, delegate to `_mapped_by_step_to_smk()` instead.

**Template:**
```python
def _task_to_rule(step, dag, rule_names):
    task = step.task or {}
    body = task.get('body', '')
    rname = _rule_name(step.id)
    has_map = step.map is not None
    lines = [f'rule {rname}:', '']

    # input section
    inputs = task.get('inputs', {})
    if inputs:
        lines.append('    input:')
        for in_name in inputs:
            binding = step.bindings.get(in_name)
            path_expr = _binding_to_path(binding, in_name, step, dag, has_map)
            lines.append(f'        {_san(in_name)}={path_expr},')
        lines.append('')

    # output section
    outputs = task.get('outputs', {})
    if outputs:
        lines.append('    output:')
        for out_name in outputs:
            default_spec = outputs[out_name].get('default')
            if default_spec:
                rendered = _interp_to_smk(default_spec)
                path_expr = repr(rendered)
            else:
                path_expr = repr(f'results/{step.id}/{out_name}')
            lines.append(f'        {_san(out_name)}={path_expr},')
        lines.append('')

    # params section (for non-file literal bindings)
    params = _collect_params(step, dag)
    if params:
        lines.append('    params:')
        for name, expr in params:
            lines.append(f'        {name}={expr},')
        lines.append('')

    # directives
    directives = _emit_resources(step)
    for d in directives:
        lines.append(f'    {d}')
    if directives:
        lines.append('')

    # shell section
    if body.strip():
        interp_body = _interpolate_shell(body, step, task)
        lines.append('    shell:')
        lines.append(f'        "{interp_body}"')
        lines.append('')

    lines.append('}')
    return '\n'.join(lines)
```

#### `_binding_to_path(binding, port_name, step, dag, has_map) → str` (~40 lines)

Convert a DAG binding value to a Snakemake file path expression.

```python
def _binding_to_path(binding, port_name, step, dag, has_map):
    if isinstance(binding, Input):
        if has_map:
            # If this step is mapped, the input varies per row → wildcard
            return repr(f'inputs/{binding.name}/{{{binding.name}}}')
        return repr(f'inputs/{binding.name}')

    if isinstance(binding, Field):
        if isinstance(binding.source, Input):
            root = binding.source
            tail = field_path_after_first(binding)
            if tail:
                return repr(f'{{inputs.{root.name}}}/{tail}')
            return f'inputs.{root.name}'
        if isinstance(binding.source, StepCall):
            prev_step = binding.source
            prev_rule = _rule_name(prev_step.id)
            out_spec = (prev_step.task or {}).get('outputs', {}).get(binding.name, {})
            default_spec = out_spec.get('default')
            if default_spec:
                rendered = _interp_to_smk(default_spec)
                return repr(rendered)
            return repr(f'results/{prev_step.id}/{binding.name}')

    if isinstance(binding, Literal):
        if isinstance(binding.value, str):
            return repr(binding.value)
        return json.dumps(binding.value)

    if isinstance(binding, Record):
        # flatten record: use a lambda input function
        fields = []
        for fname in sorted(binding.fields.keys()):
            fpath = _binding_to_path(binding.fields[fname], fname, step, dag, has_map)
            fields.append(f'{fname}: {fpath}')
        return '{{' + ', '.join(fields) + '}}'

    raise ValueError(f'Unsupported binding for Snakemake: {type(binding).__name__}')
```

#### `_san(name) → str` (~3 lines)

Sanitize a port name for use as a Snakemake named-input key.

```python
def _san(name):
    return name.replace('-', '_')
```

#### `_interp_to_smk(value) → str` (~10 lines)

Render a SWL `word` interpolation dict into a Snakemake format string using `word_interp()` from `common.py`.

```python
def _interp_to_smk(value):
    return word_interp(
        value,
        literal_fn=lambda text: text,
        var_fn=lambda name: f'{{{name}}}',
        expr_fn=lambda text: f'${{{{{text}}}}}',
    )
```

- `literal` text → passed through verbatim.
- `var` (e.g. `sample_id`) → `{sample_id}` (Snakemake wildcard / python format).
- `expr` (e.g. `"out_" + sample`) → `${{"out_" + sample}}` (double-braced escaped Python expression).

#### `_interpolate_shell(body, step, task) → str` (~20 lines)

Rewrite bash shell body to use Snakemake `{input.port}`, `{output.port}`, `{params.port}` references.

```python
def _interpolate_shell(body, step, task):
    input_names = set(task.get('inputs', {}).keys())
    output_names = set(step.outputs)
    param_names = set(p[0] for p in _collect_params(step, None))

    result = []
    for line in body.split('\n'):
        # Replace $var or ${var} references
        # (simple heuristic: match $name for known ports)
        replaced = _replace_shell_vars(line, input_names, output_names, param_names)
        result.append(replaced)
    return '\n'.join(result)
```

The replacer substitutes `$fastq1` → `{input.fastq1}`, `$bam` → `{output.bam}`, and `$threshold` → `{params.threshold}` for known names.

#### `_collect_params(step, dag) → list` (~15 lines)

Collect non-file literal bindings into `params:` entries:

```python
def _collect_params(step, dag):
    params = []
    task_inputs = (step.task or {}).get('inputs', {})
    for in_name, spec in task_inputs.items():
        typ = spec.get('type', 'str')
        if typ in ('str', 'int', 'float', 'memory', 'time'):
            binding = step.bindings.get(in_name)
            if isinstance(binding, Literal):
                params.append((_san(in_name), json.dumps(binding.value)))
            elif isinstance(binding, Input):
                params.append((_san(in_name), f'config[{json.dumps(binding.name)}]'))
    return params
```

#### `_emit_resources(step) → list[str]` (~20 lines)

Map SWL `run` parameters to Snakemake directives using `run_value()` from `common.py`:

```python
def _emit_resources(step):
    run = (step.task or {}).get('run', {})
    directives = []

    cpu = run_value(run, 'cpu')
    if cpu is not None:
        directives.append(f'threads: {cpu}')

    memory = run_value(run, 'memory')
    if memory is not None:
        directives.append(f'resources:')
        directives.append(f'    mem_mb={memory}')

    time_val = run_value(run, 'time')
    if time_val is not None:
        directives.append(f'resources:')
        directives.append(f'    runtime_minutes={time_val}')

    image = run_value(run, 'image')
    if image is not None:
        directives.append(f'container: "docker://{image}"')

    return directives
```

#### `_dag_to_smk(dag, workflow_id, rule_names) → str` (~30 lines)

Generate the workflow body: a `rule all:` with final outputs as targets. For non-mapped (simple) workflows, each step call is wired via filename matching. The `rule all:` is the entry point that Snakemake uses to determine what files to produce.

```python
def _dag_to_smk(dag, workflow_id, rule_names):
    lines = [f'# Workflow: {_wf_name(workflow_id)}', '']

    # Wrap inputs into a config reference
    if dag.inputs:
        lines.append('# Workflow inputs — define in config.yaml or via --config')
        lines.append('')

    # rule all: — target rule
    lines.append('rule all:')
    lines.append('    input:')
    for name, output in dag.outputs.items():
        value = output.value if isinstance(output, OutputSpec) else output
        path_expr = _output_to_path(name, value, dag)
        lines.append(f'        {_san(name)}={path_expr},')
    lines.append('')

    return '\n'.join(lines)
```

#### `_output_to_path(name, value, dag) → str` (~15 lines)

Convert a workflow output value to a Snakemake-compatible file path:

```python
def _output_to_path(name, value, dag):
    if isinstance(value, Input):
        return repr(f'inputs/{value.name}')
    if isinstance(value, Literal):
        return repr(value.value)
    if isinstance(value, Field):
        if isinstance(value.source, StepCall):
            prev = value.source
            spec = (prev.task or {}).get('outputs', {}).get(value.name, {})
            default = spec.get('default')
            if default:
                rendered = _interp_to_smk(default)
                return repr(rendered)
            return repr(f'results/{prev.id}/{value.name}')
        if isinstance(value.source, Input):
            return repr(f'{{inputs.{value.source.name}}}/{value.name}')
        raise ValueError(f'Unsupported Field source for output: {type(value.source).__name__}')
    if isinstance(value, Record):
        # Decompose record output into multiple rule all inputs
        parts = []
        for fname, fval in value.fields.items():
            parts.append(_output_to_path(f'{name}_{fname}', fval, dag))
        return ', '.join(parts)
    raise ValueError(f'Unsupported output value for Snakemake: {type(value).__name__}')
```

#### `_mapped_step_to_rule(step, dag, rule_names) → str` (~35 lines)

For SWL `map f xs` steps: generate a rule with wildcard-based scatter.

The key insight: Snakemake's wildcard system maps directly to SWL's `map`. The table columns become wildcards in the filenames.

```python
def _mapped_step_to_rule(step, dag, rule_names):
    map_info = step.map or {}
    source = map_info.get('source', {})
    schema = step.input_schema or {}
    task = step.task or {}

    # Determine wildcard names from scatter ports or schema
    scatter_ports = map_info.get('scatter', [])
    wildcard_names = scatter_ports if scatter_ports else list(schema.keys())

    # Build output paths with wildcards
    output_specs = task.get('outputs', {})
    for out_name in step.outputs:
        spec = output_specs.get(out_name, {})
        default_spec = spec.get('default')

    # Build input paths with wildcard substitution
    # For table-sourced columns, column_input_name() maps column → input name
    columns = table_columns(source)

    ...
```

Full implementation follows the same structure as `_task_to_rule()` but adds wildcards to output filenames and uses `{wildcards.name}` references in input path expressions.

#### `_mapped_by_step_to_rule(step, dag, rule_names) → str` (~80 lines)

For SWL `map_by f key xs`: generate two rules — a grouping rule and a processing rule.

**Grouping rule:** Aggregates data per key using Snakemake's `group_by()` function (7.x+) or a Python `run:` block.

```python
rule group_{step_id}:
    input:
        # All input files by wildcard
        ...
    output:
        groups="{key}/data.txt"
    run:
        # Python: read input, group by key, write per-group files
        ...
```

**Processing rule:** Uses the grouping key as a wildcard:

```python
rule {step_id}:
    input:
        data="{{key}}/data.txt"
    output:
        result="{{key}}/result.txt"
    shell:
        "tool {input.data} > {output.result}"
```

The implementation uses `source_kind()`, `source_input_name()`, and `table_columns()` from `common.py` to determine table structure, then emits appropriate `expand()` and `group_by()` calls.

#### `_subworkflow_to_smk(step, parent_id) → str` (~20 lines)

Recursively transpile sub-workflow steps into a module Snakefile. Uses the same `transpile_dag_dict()` entry point (similar to `transpile/wdl/emit.py:_subworkflow_to_wdl`).

```python
def _subworkflow_to_smk(step, parent_id):
    definition = step.task or {}
    dag_data = definition.get('dag', {})
    if not dag_data:
        return ''
    wf_id = f'{parent_id}_{step.id}'
    return transpile_dag_dict(dag_data, workflow_id=wf_id, _top_level=False)
```

At the call site in the parent workflow, use Snakemake's `module` directive:

```python
module {step.id}:
    snakefile:
        "{step.id}.smk"

use rule * from {step.id} as {step.id}_*
```

Or inline all rules with prefixed names (simpler, matching the pattern in `transpile/nf/emit.py`).

### Estimated effort by component

| Component | Lines | Uses from `common.py` | Difficulty |
|-----------|-------|----------------------|-----------|
| `transpile_smk_dict()` (main entry point) | 30 | — | Easy |
| `_rule_name()` / `_wf_name()` / `_san()` | 6 | `step_name`, `workflow_name`, `emit_name` | Trivial |
| `_task_to_rule()` — basic step → rule | 80 | `run_value` | Easy |
| `_binding_to_path()` — value-to-path expression | 40 | `field_chain_parts`, `field_path_after_first` | Medium |
| `_interp_to_smk()` — string interpolation | 10 | `word_interp` | Easy |
| `_interpolate_shell()` — bash rewriting | 20 | — | Easy |
| `_collect_params()` — params section | 15 | — | Easy |
| `_emit_resources()` — threads, resources, container | 20 | `run_value` | Easy |
| `_dag_to_smk()` — rule all + wiring | 30 | — | Medium |
| `_output_to_path()` — output value → file path | 15 | `field_chain_parts` | Medium |
| `_mapped_step_to_rule()` — scatter via wildcards | 35 | `source_kind`, `source_input_name`, `table_columns`, `column_input_name` | Easy |
| `_mapped_by_step_to_rule()` — grouped scatter | 80 | `source_kind`, `source_input_name`, `table_columns`, `column_input_name` | Medium |
| `_subworkflow_to_smk()` — module/inline | 20 | — | Medium |
| `transpile_dag_file()` — file entry point | 3 | — | Trivial |
| `cli.py`, `__init__.py`, `__main__.py` | 12 | `swl.transpile._cli.run` | Trivial |
| **Total** | **~416** | | |

### Unsupported or complex cases

1. **Literal workflow outputs**: Snakemake does not support constant/literal outputs at the workflow level (same limitation as WDL and CWL). Raise a clear error, matching existing transpiler behavior.

2. **Merge bindings in outputs**: Like WDL, Snakemake cannot express record merge in workflow outputs. These must be flattened by the compiler's finalization phase before transpilation.

3. **Expression-based output paths**: SWL output paths with `${expr}` interpolation need a `params` or `run:` Python block to compute the path at runtime. Feasible but adds complexity.

4. **Non-file table sources**: If a `map` source is in-memory (not filesystem files), the transpiler must generate a CSV/JSON file and use Snakemake `params` to pass values — similar to CWL's approach with `ExpressionTool`.

---

## Comparison with Existing Transpilers

| Aspect | CWL | WDL | Nextflow | Snakemake (estimated) |
|--------|-----|-----|----------|-----------------------|
| Lines of code | 572 | 520 | 344 | ~520 |
| map_by support | JS grouping tool | `collect_by_key` | `groupTuple` | `group_by()` or Python |
| Record bindings | `ExpressionTool` | `struct` | `tuple` | Flatten to named inputs |
| Sub-workflows | Nested `$graph` | Nested `workflow` | `workflow` module | `module` directive |
| Dependency model | Port-based | Port-based | Channel-based | Filename-based |
| Difficulty of mapping record/merge | Hard | Hard | Medium | Medium |

Snakemake's filename-based model is the biggest conceptual gap from SWL's port-based model. However, using consistently generated output filename patterns across rules makes this tractable. The Nextflow transpiler faces a similar channel-to-filename bridging challenge and succeeds in 344 lines, suggesting the Snakemake transpiler should be achievable in ~520 lines.

---

## Risk Assessment

- **Low risk**: Simple rules, resource requirements, string interpolation, named ports, scatter (wildcards)
- **Medium risk**: `map_by` (grouped scatter), sub-workflow modularization, record bindings
- **High risk**: None identified. All SWL constructs have a reasonable Snakemake counterpart.

---

## Conclusion

A Snakemake transpiler is **feasible and implementable** with approximately 500 lines of Python, following the existing transpiler pattern in `transpile/cwl/`, `transpile/wdl/`, and `transpile/nf/`. The main effort will be in generating consistent filename patterns for Snakemake's wildcard-based dependency resolution and handling `map_by` grouped scatter. No architectural changes to the compiler's DAG output format are needed — the transpiler operates on the existing normalized DAG contract.
