# WDL Transpiler: Implementation Plan

Transpile SWL compiled DAG JSON to WDL 1.1. Follows the same pattern as `swl/transpile/cwl/` — reads DAG JSON, emits target code, never reads source files.

---

## 1. WDL 1.1 Concepts

| Concept | Role | SWL/DAG equivalent |
|---------|------|--------------------|
| `task` | A bash script encapsulated as a reusable function | `StepCall` with `task.body` |
| `workflow` | Composition block wiring task/sub-workflow calls | `DAG.inputs` + `DAG.steps` + `DAG.outputs` |
| `input {}` | Declares task or workflow inputs | `InputSpec` / `ParamSpec` |
| `output {}` | Declares task or workflow outputs | `ParamSpec` with glob/default |
| `command <<< >>>` | Bash script with `~{}` placeholders | `step.task.body` |
| `call task { input: x = y }` | Invoke a task with wired inputs | `StepCall.bindings` |
| `scatter (x in arr) { ... }` | Parallel execution over array elements | `MappedStep` (simple `map`) |
| `requirements {}` | Resource requirements (container, cpu, memory) | `RunParamSpec` |
| `struct` | User-defined record type | `Record` binding shape |
| `collect_by_key()` | Group array elements by key column | `map_by` grouping |
| `after` | Explicit dependency ordering | Step dependency edges |
| `Array[X]` | Array type | `[file]`, `[str]`, `[int]`, `[float]` |
| `Pair[X,Y]` | Tuple type | Record merge internals |
| `None` | Optional value marker | Unset optional parameter |
| `import` | Include another WDL file | Sub-workflow `definition.dag` |

---

## 2. Transpiler Structure

```
swl/transpile/wdl/
    __init__.py        # package marker
    cli.py             # CLI entry point: wdl transpile <dag.json> -o <output.wdl>
    __main__.py         # python -m swl.transpile.wdl
    emit.py             # core transpilation logic
    test_emit.py        # tests
```

### Entry point

```python
# cli.py
import argparse, json
from swl.ir.dag import DAG
from swl.transpile.wdl.emit import transpile_dag_dict

def main():
    ap = argparse.ArgumentParser('Transpile compiled DAG JSON to WDL 1.1')
    ap.add_argument('input', help='compiled DAG json path')
    ap.add_argument('-o', '--output', help='output .wdl path')
    args = ap.parse_args()
    data = json.load(open(args.input))
    wdl = transpile_dag_dict(data)
    if args.output:
        with open(args.output, 'w') as f:
            f.write(wdl)
    else:
        print(wdl)
```

### `emit.py` structure

```python
def transpile_dag_dict(data, workflow_id='main'):
    dag = DAG.from_dict(data)
    _validate_supported(dag)

    structs = _collect_structs(dag)   # gather record shapes for struct definitions
    tasks = {}
    for step in dag.steps:
        if step.id not in tasks:
            if step.type == 'workflow':
                tasks[step.id] = _subworkflow_to_wdl(step, workflow_id)
            else:
                tasks[step.id] = _task_to_wdl(step)

    workflow = _dag_to_wdl(dag, workflow_id, tasks)

    lines = []
    lines.append('version 1.1\n')
    # Emit struct definitions
    for s in structs:
        lines.append(s)
        lines.append('')
    # Emit task definitions
    for name, body in tasks.items():
        lines.append(body)
        lines.append('')
    # Emit workflow block
    lines.append(workflow)
    return '\n'.join(lines)
```

---

## 3. DAG-to-WDL Mapping

### 3.1 Task step → WDL task

A `StepCall` with `type: "task"` becomes a WDL `task` block.

```wdl
task align {
    input {
        File fastq1
        File fastq2
        File ref
        String outbase
    }

    command <<<
        bwa mem -t ~{cpu} ~{ref} ~{fastq1} ~{fastq2} |
            samtools view -b - > ~{outbase}.bam
    >>>

    output {
        File bam = "~{outbase}.bam"
    }

    requirements {
        container: "djhshih/seqkit:0.1"
        cpu: 2
        memory: "8192 MB"
    }
}
```

**Mapping rules:**

| DAG field | WDL 1.1 |
|-----------|---------|
| `id` | Task name (sanitized: lowercase, underscores) |
| `run.cpu.value` | `cpu: <value>` |
| `run.memory.value` | `memory: '<value> MB'` |
| `run.image.value` | `container: '<value>'` |
| `run.time.value` | `time_minutes: <value>` |
| `inputs` | `input {}` block, each param |
| `outputs` | `output {}` block, each param |
| `outputs[].default` | File path in `File` declaration with interpolation |
| `task.body` | `command <<< >>>` block |
| `bindings` | Expression wiring in `workflow` block's `call { input: ... }` |

**Type mapping:**

| SWL type | WDL type |
|----------|----------|
| `file` | `File` |
| `str` | `String` |
| `int` | `Int` |
| `float` | `Float` |
| `[file]` | `Array[File]` |
| `[str]` | `Array[String]` |
| `[int]` | `Array[Int]` |
| `[float]` | `Array[Float]` |

```python
def _wdl_type(swl_type, optional=False):
    base = {
        'file': 'File',
        'str': 'String',
        'int': 'Int',
        'float': 'Float',
        '[file]': 'Array[File]',
        '[str]': 'Array[String]',
        '[int]': 'Array[Int]',
        '[float]': 'Array[Float]',
    }.get(swl_type, 'String')
    if optional:
        return base + '?'
    return base
```

**Task name sanitization:**

```python
def _task_name(step_id):
    name = step_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'task'
    if name[0].isdigit():
        name = '_' + name
    return name.lower()
```

### 3.2 Workflow block

The top-level workflow block contains `input {}`, `call` statements, scatters, and `output {}`.

```wdl
workflow main {
    input {
        File fastq1
        File fastq2
        File ref
        String outbase
    }

    call align {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            ref = ref,
            outbase = outbase
    }

    output {
        File bam = align.bam
    }
}
```

**Workflow input declaration:**

| Input type | WDL declaration |
|------------|-----------------|
| `file` | `File <name>` |
| `str` | `String <name>` |
| `int` | `Int <name>` |
| `file` array | `Array[File] <name>` |
| optional `file` | `File? <name>` |

**Binding-to-call-expression wiring:**

Each binding in `step.bindings` maps to a named input in the `call` block:

| Binding `source` | WDL call input expression |
|---|---|
| `input` | `<input_name>` |
| `step_output` / `step` | `<STEP_NAME>.<output_name>` |
| `literal` | literal value |
| `field` | `<expr>.<field_name>` |
| `merge` | Flatten to individual fields (see §10) |
| `record` | `StructName { fields... }` (struct literal) |

```python
def _dag_to_wdl(dag, workflow_id, tasks):
    lines = [f'workflow {_wf_name(workflow_id)} {{', '']

    # Input section
    if dag.inputs:
        lines.append('    input {')
        for name, spec in dag.inputs.items():
            lines.append(f'        {_wdl_type(spec.type)} {name}')
        lines.append('    }')
        lines.append('')

    # Step calls in topological order
    for step in dag.steps:
        if getattr(step, 'map', None) is not None:
            lines.extend(_mapped_step_to_wdl(step, tasks))
            continue

        tname = _task_name(step.id)
        alias = _call_alias(step.id)
        if alias == tname:
            lines.append(f'    call {tname} {{')
        else:
            lines.append(f'    call {tname} as {alias} {{')
        lines.append('        input:')
        for in_name, binding in step.bindings.items():
            expr = _binding_to_wdl_expr(binding, step.id, dag)
            lines.append(f'            {in_name} = {expr},')
        lines.append('    }')
        lines.append('')

    # Output section
    if dag.outputs:
        lines.append('    output {')
        for name, binding in dag.outputs.items():
            expr = _binding_to_wdl_expr(binding, None, dag)
            typ = _infer_output_type(name, binding, dag)
            lines.append(f'        {typ} {name} = {expr}')
        lines.append('    }')

    lines.append('}')
    return '\n'.join(lines)
```

### 3.3 Binding-to-expression

```python
def _binding_to_wdl_expr(binding, current_step_id, dag):
    if isinstance(binding, Input):
        return binding.name

    if isinstance(binding, Literal):
        return _literal_to_wdl(binding.value)

    if isinstance(binding, Field):
        source_expr = _binding_to_wdl_expr(binding.source, current_step_id, dag)
        return f'{source_expr}.{binding.name}'

    if isinstance(binding, Merge):
        raise ValueError(
            f'Merge bindings must be flattened before WDL transpilation'
        )

    if isinstance(binding, Record):
        struct_name = _struct_name_for(binding)
        fields = []
        for fname, fbinding in binding.fields.items():
            fexpr = _binding_to_wdl_expr(fbinding, current_step_id, dag)
            fields.append(f'{fname}: {fexpr}')
        return f'{struct_name} {{{", ".join(fields)}}}'

    if isinstance(binding, (StepCall, MappedStep)):
        return _call_alias(binding.id)

    raise ValueError(f'Unsupported binding for WDL: {type(binding).__name__}')


def _literal_to_wdl(value):
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if value is None:
        return 'None'
    return str(value)
```

---

## 4. Mapped Steps (`map`) → Scatter

A `MappedStep` without `group_by` represents `map f xs`. In WDL this maps to a `scatter` block.

### Strategy: scatter over zipped array columns

The mapped step's `input_schema` defines per-row types. The `map.source` identifies the table source. Columns of the table are `Array[X]`-typed workflow inputs.

```wdl
workflow main {
    input {
        Array[File] fastq1
        Array[File] fastq2
        String outbase      # broadcast: same value for all rows
    }

    scatter (idx in range(length(fastq1))) {
        call align {
            input:
                fastq1 = fastq1[idx],
                fastq2 = fastq2[idx],
                outbase = outbase
        }
    }

    output {
        Array[File] bam = align.bam
    }
}
```

**Key decisions:**
1. Scatter over indices (`range(length(...))`) rather than zipping directly, because broadcast inputs (single values shared across all rows) must be handled without creating per-row arrays.
2. The scatter body can contain only `call` and `if` statements. Intermediate expressions must be embedded in call inputs or declared before the scatter.
3. WDL automatically gathers scattered outputs into `Array[]` — no explicit gather step needed.

**Identifying scatter inputs vs broadcast inputs:**
- Inputs that come from `map.source` (table columns) are scatter inputs → indexed with `[idx]`
- Inputs that come from elsewhere (literals, step outputs from non-mapped steps) are broadcast → used directly

```python
def _mapped_step_to_wdl(step, tasks):
    lines = []
    tname = _task_name(step.id)
    alias = _call_alias(step.id)
    s_var = f'{step.id}_i'

    source = step.map.get('source', {})
    table_columns = _get_table_columns(source)
    len_expr = _derive_length_expr(source, step.bindings, table_columns)
    lines.append(f'    scatter ({s_var} in range({len_expr})) {{')

    call_line = f'        call {tname}'
    if alias != tname:
        call_line += f' as {alias}'
    lines.append(call_line + ' {')
    lines.append('            input:')
    for in_name in step.task.get('inputs', {}).keys():
        binding = step.bindings.get(in_name)
        if binding is None:
            if in_name in table_columns:
                expr = f'{_column_input_name(source, in_name)}[{s_var}]'
            else:
                expr = in_name
        elif isinstance(binding, Input):
            expr = f'{binding.name}[{s_var}]' if binding.name in table_columns else binding.name
        elif isinstance(binding, Literal):
            expr = _literal_to_wdl(binding.value)
        elif isinstance(binding, Field):
            base = _binding_to_wdl_expr(binding.source, step.id, None)
            if isinstance(binding.source, (StepCall, MappedStep)):
                expr = f'{base}.{binding.name}[{s_var}]'
            else:
                expr = f'{base}.{binding.name}'
        else:
            expr = _binding_to_wdl_expr(binding, step.id, None)
        lines.append(f'                {in_name} = {expr},')
    lines.append('        }')
    lines.append('    }')
    lines.append('')

    return lines


def _get_table_columns(source):
    if source.get('source') == 'table':
        return set(source.get('columns', {}).keys())
    return set()


def _column_input_name(source, col_name):
    if source.get('source') == 'input':
        return source['name']
    if source.get('source') == 'table':
        col = source.get('columns', {}).get(col_name, {})
        if isinstance(col, dict) and col.get('source') == 'input':
            return col['name']
    return col_name


def _derive_length_expr(source, bindings, table_columns):
    if source.get('source') == 'input':
        return f'length({source["name"]})'
    if source.get('source') == 'table':
        columns = source.get('columns', {})
        col_names = list(columns.keys())
        if col_names:
            col = columns[col_names[0]]
            if isinstance(col, Input):
                return f'length({col.name})'
            if isinstance(col, Field) and isinstance(col.source, (StepCall, MappedStep)):
                return f'length({_call_alias(col.source.id)}.{col.name})'
    for name, binding in bindings.items():
        if isinstance(binding, Input):
            return f'length({binding.name})'
    return '1'
```

---

## 5. `map_by` → `collect_by_key()` + Scatter

WDL 1.1 has `collect_by_key()` — a standard library function that groups array elements by a key column. SWL's `map_by` maps directly to this.

### Strategy: `collect_by_key()` with per-group scatter

For `map_by align "sample" samples` where `samples` has columns `sample, fastq1, fastq2`:

```wdl
workflow main {
    input {
        Array[String] sample
        Array[File] fastq1
        Array[File] fastq2
    }

    # 1. Zip columns into pairs, collect by key
    Array[Pair[String, Pair[File, File]]] grouped = collect_by_key(
        zip(sample, zip(fastq1, fastq2))
    )

    # 2. Scatter over groups; each group element is a Pair(key, Array[values])
    scatter (group in grouped) {
        String group_sample = group.left
        Array[File] group_fastq1 = group.right.left
        Array[File] group_fastq2 = group.right.right

        call align_group {
            input:
                sample = group_sample,
                fastq1 = group_fastq1,
                fastq2 = group_fastq2
        }
    }
}
```

**Simplified approach for multi-column tables:**

For tables with many columns, use an intermediate grouping step (inline declarations + `collect_by_key`) before the scatter:

```wdl
    # Zip all table columns
    Array[Pair[String, Pair[File, Pair[File, String]]]] zipped = zip(
        sample,
        zip(fastq1, zip(fastq2, outbase))
    )

    # Group by sample (first element)
    Array[Pair[String, Array[Pair[File, Pair[File, String]]]]] grouped = collect_by_key(zipped)

    scatter (g in grouped) {
        String key = g.left
        Array[File] g_fastq1 = select_all(map(g.right, r -> r.left))
        Array[File] g_fastq2 = select_all(map(g.right, r -> r.right.left))
        Array[String] g_outbase = select_all(map(g.right, r -> r.right.right))

        call align_group {
            input:
                fastq1 = g_fastq1,
                fastq2 = g_fastq2,
                outbase = g_outbase
        }
    }
```

### Emitting `map_by`

```python
def _mapped_by_step_to_wdl(step, tasks):
    lines = []
    tname = _task_name(step.id)
    alias = _call_alias(step.id)
    map_info = step.map
    group_key = map_info.get('group_by')
    source = map_info.get('source', {})
    input_schema = step.input_schema or {}

    # Collect column binding names in order
    col_names = list(input_schema.keys())
    if group_key and group_key not in col_names:
        col_names = [group_key] + col_names

    # Build zip expression for all columns
    # zip(a, zip(b, zip(c, d)))  -- repeatedly pair
    zip_expr = _build_zip_chain(col_names, source)
    grouped_var = f'{step.id}_grouped'
    lines.append(f'    Array[Pair[{_key_type(group_key, input_schema)}, Array[{_val_type(col_names, input_schema)}]]] {grouped_var} = collect_by_key({zip_expr})')
    lines.append('')

    # Scatter over groups
    g_var = f'{step.id}_g'
    lines.append(f'    scatter ({g_var} in {grouped_var}) {{')

    # Extract key and value arrays inside scatter
    lines.append(f'        {_wdl_type(input_schema.get(group_key, "str"))} {group_key}_val = {g_var}.left')
    for col in col_names:
        if col == group_key:
            continue
        t = _wdl_type(input_schema.get(col, 'str'))
        lines.append(f'        Array[{t}] {col}_vals = {g_var}.right.{_col_access_path(col, col_names)}')

    # Call the grouped process
    call_line = f'        call {tname}'
    if alias != tname:
        call_line += f' as {alias}'
    lines.append(call_line + ' {')
    lines.append('            input:')
    for in_name in step.task.get('inputs', {}).keys():
        if in_name == group_key:
            lines.append(f'                {in_name} = {group_key}_val,')
        elif in_name in col_names:
            lines.append(f'                {in_name} = {in_name}_vals,')
        else:
            binding = step.bindings.get(in_name)
            if binding:
                expr = _binding_to_wdl_expr(binding, step.id, None)
                lines.append(f'                {in_name} = {expr},')
            else:
                lines.append(f'                {in_name} = {in_name},')
    lines.append('        }')
    lines.append('    }')
    lines.append('')

    return lines


def _build_zip_chain(col_names, source):
    if not col_names:
        return '[]'
    inner = _column_expr(col_names[-1], source)
    for col in reversed(col_names[:-1]):
        col_expr = _column_expr(col, source)
        inner = f'zip({col_expr}, {inner})'
    return inner


def _column_expr(col, source):
    if source.get('source') == 'input':
        return source['name']
    columns = source.get('columns', {})
    col_binding = columns.get(col, {})
    if isinstance(col_binding, dict) and col_binding.get('source') == 'input':
        return col_binding['name']
    return col


def _key_type(key, schema):
    return _wdl_type(schema.get(key, 'str'))


def _val_type(col_names, schema):
    if len(col_names) == 1:
        return _wdl_type(schema.get(col_names[0], 'str'))
    if len(col_names) == 2:
        return f'Pair[{_wdl_type(schema.get(col_names[0], "str"))}, {_wdl_type(schema.get(col_names[1], "str"))}]'
    inner = _wdl_type(schema.get(col_names[-1], 'str'))
    for col in reversed(col_names[1:-1]):
        inner = f'Pair[{_wdl_type(schema.get(col, "str"))}, {inner}]'
    return f'Pair[{_wdl_type(schema.get(col_names[0], "str"))}, {inner}]'


def _col_access_path(col, col_names):
    """Return the .left/.right path to extract a non-key column from the Pair chain."""
    idx = col_names.index(col)
    if idx == 0:
        return 'left'
    path = 'right'
    for _ in range(1, idx):
        path += '.right'
    path += '.left'
    return path
```

---

## 6. Sub-workflow Steps

A step with `type: "workflow"` has a `definition` containing a recursive DAG. WDL 1.1 supports calling sub-workflows.

**Strategy:** Emit all sub-workflows in the same file as additional `workflow` blocks. WDL 1.1 allows multiple workflow blocks; only the first is the entry point.

```wdl
# Parent workflow:
call sub_pipeline as sub {
    input:
        fastq1 = fastq1,
        fastq2 = fastq2
}
```

**Recursive emission:**

```python
def _subworkflow_to_wdl(step, parent_id):
    definition = step.task or {}
    dag_data = definition.get('dag', {})
    if not dag_data:
        return ''
    wf_id = f'{parent_id}_{step.id}'
    return transpile_dag_dict(dag_data, workflow_id=wf_id)
```

The parent workflow's call references the sub-workflow by the `_wf_name` of the recursive emission. Since all workflow blocks are in the same file, WDL engines resolve them by name.

```python
def _wf_name(workflow_id):
    name = workflow_id.replace('-', '_').lstrip('_')
    if not name:
        name = 'main'
    if name[0].isdigit():
        name = '_' + name
    return name
```

---

## 7. Resource and Runtime Mapping

WDL 1.1 uses the `requirements {}` section for runtime resources. The `runtime {}` section is also still valid but `requirements` is the canonical location.

| DAG field | WDL 1.1 requirement |
|-----------|---------------------|
| `run.cpu.value` | `cpu: <value>` |
| `run.memory.value` | `memory: '<value> MB'` |
| `run.image.value` | `container: '<value>'` |
| `run.time.value` | `time_minutes: <value>` |

```python
def _emit_requirements(step):
    attrs = []
    run = step.run if isinstance(step, (StepCall, MappedStep)) else {}

    for name, spec in run.items():
        if isinstance(spec, dict):
            value = spec.get('value')
        elif hasattr(spec, 'value'):
            value = spec.value
        else:
            continue
        if value is None:
            continue

        if name == 'cpu':
            attrs.append(f'        cpu: {value}')
        elif name == 'memory':
            attrs.append(f'        memory: "{value} MB"')
        elif name == 'image':
            attrs.append(f'        container: "{value}"')
        elif name == 'time':
            attrs.append(f'        time_minutes: {value}')

    if not attrs:
        return ''
    return '    requirements {\n' + '\n'.join(attrs) + '\n    }'
```

---

## 8. Interpolation Handling

WDL 1.1 uses `~{}` syntax for expression placeholders in `command <<< >>>` blocks.

```python
def _interp_to_wdl(value):
    if value is None:
        return None
    if value.get('kind') == 'word':
        parts = value.get('parts', [])
        result = ''
        for part in parts:
            if part.get('kind') == 'literal':
                result += part['text']
            elif part.get('kind') == 'var':
                result += f"~{{{part['name']}}}"
            elif part.get('kind') == 'expr':
                result += f"~{{{part['text']}}}"
        return result
    return None
```

**Usage in output declarations:**

```python
def _emit_outputs(step):
    lines = []
    lines.append('    output {')
    for name, spec in step.task.get('outputs', {}).items():
        typ = _wdl_type(spec.get('type'))
        default = spec.get('default')
        if default:
            path_expr = _interp_to_wdl(default)
            lines.append(f'        {typ} {name} = "{path_expr}"')
        else:
            lines.append(f'        {typ} {name} = "*.{name}"')
    lines.append('    }')
    return '\n'.join(lines)
```

---

## 9. Workflow Output Wiring

Output bindings at the `DAG.outputs` level are converted to `output {}` declarations in the workflow block.

| Output binding | WDL expression |
|---|---|
| `{step: "align", output: "bam"}` | `align.bam` |
| `{source: "input", name: "fastq1"}` | `fastq1` |
| `{source: "literal", value: 42}` | `42` |
| `{source: "field", field: "x", value: ...}` | `<expr>.x` |

When the output of a mapped step is referenced, WDL's scatter-gather semantics automatically collect it into an `Array`.

```python
def _infer_output_type(name, binding, dag):
    if isinstance(binding, Field) and isinstance(binding.source, (StepCall, MappedStep)):
        step = binding.source
        out_type = step.task.get('outputs', {}).get(binding.name, {}).get('type', 'str')
        wdl_type = _wdl_type(out_type)
        if getattr(step, 'map', None) is not None:
            return f'Array[{wdl_type}]'
        return wdl_type
    if isinstance(binding, Input):
        spec = dag.inputs.get(binding.name)
        return _wdl_type(spec.type if spec else 'str')
    if isinstance(binding, Literal):
        return _infer_literal_type(binding.value)
    return 'String'
```

---

## 10. Merge and Record Handling

### Record Merge

WDL has no record merge operator. Merge bindings (`{source: "merge", ...}`) must be **flattened before reaching the transpiler**.

```python
def _validate_binding(value, step_id):
    if isinstance(value, Merge):
        raise ValueError(
            f'WDL transpilation does not support merge bindings ({step_id}). '
            f'Merge bindings must be flattened before transpilation.'
        )
```

The flattening should be done by the compiler (in `force.py`) before DAG output. Currently, chain desugaring (`A | B | C`) already flattens merges in most cases.

### Record Binding → Struct

Record bindings (`{source: "record", fields: {...}}`) map to WDL 1.1 `struct` literals. A `struct` definition must be emitted for each unique record shape.

```wdl
struct AlignInput {
    File fastq1
    File fastq2
    String outbase
}

# Usage:
call align {
    input:
        params = AlignInput {
            fastq1: fastq1_ch,
            fastq2: fastq2_ch,
            outbase: "output"
        }
}
```

**Struct collection:**

```python
def _collect_structs(dag):
    """Walk all bindings and collect unique record shapes for struct generation."""
    structs = []
    seen = set()
    for step in dag.steps:
        for binding in step.bindings.values():
            if isinstance(binding, Record):
                shape = _record_shape(binding)
                if shape not in seen:
                    seen.add(shape)
                    structs.append(_emit_struct(shape))
    return structs


def _record_shape(binding):
    """Produce a canonical name from the field set."""
    if not isinstance(binding, Record):
        return None
    keys = tuple(sorted(binding.fields.keys()))
    return keys


def _struct_name_for(binding):
    keys = tuple(sorted(binding.fields.keys()))
    return '_Rec_' + '_'.join(k.capitalize() for k in keys)


def _emit_struct(shape):
    """Emit a struct definition string."""
    lines = [f'struct {_struct_name_for(shape)} {{']
    for name in shape:
        lines.append(f'    String {name}')
    lines.append('}')
    return '\n'.join(lines)
```

---

## 11. Complete Example

### Input DAG (`align_merge.swl` compiled)

```json
{
  "version": "1.0",
  "inputs": {
    "fastq1": {"type": "[file]", "desc": "read 1"},
    "fastq2": {"type": "[file]", "desc": "read 2"},
    "ref": {"type": "file", "desc": "reference"},
    "outbase": {"type": "str", "desc": "output base"}
  },
  "steps": [
    {
      "id": "align",
      "type": "mapped_task",
      "map": {"source": {"source": "input", "name": "fastq1"}, "group_by": null},
      "input_schema": {"fastq1": "file", "fastq2": "file", "ref": "file", "outbase": "str"},
      "output_schema": {"bam": "file"},
      "deps": [],
      "inputs": {"fastq1": {"type": "file"}, "fastq2": {"type": "file"}, "ref": {"type": "file"}, "outbase": {"type": "str"}},
      "bindings": {"ref": {"source": "input", "name": "ref"}, "outbase": {"source": "input", "name": "outbase"}},
      "outputs": {"bam": {"type": "file", "default": {"kind": "word", "parts": [{"kind": "literal", "text": ""}, {"kind": "var", "name": "outbase"}, {"kind": "literal", "text": ".bam"}]}}},
      "run": {"cpu": {"type": "int", "value": 2}, "memory": {"type": "memory", "value": 8192}, "image": {"type": "str", "value": "djhshih/seqkit:0.1"}},
      "script": "bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} > ${outbase}.bam"
    }
  ],
  "outputs": {"bam": {"step": "align", "output": "bam"}}
}
```

### Output WDL 1.1

```wdl
version 1.1

task align {
    input {
        File fastq1
        File fastq2
        File ref
        String outbase
    }

    command <<<
        bwa mem -t ~{cpu} ~{ref} ~{fastq1} ~{fastq2} > ~{outbase}.bam
    >>>

    output {
        File bam = "~{outbase}.bam"
    }

    requirements {
        cpu: 2
        memory: "8192 MB"
        container: "djhshih/seqkit:0.1"
    }
}

workflow main {
    input {
        Array[File] fastq1
        Array[File] fastq2
        File ref
        String outbase
    }

    scatter (align_i in range(length(fastq1))) {
        call align as align {
            input:
                fastq1 = fastq1[align_i],
                fastq2 = fastq2[align_i],
                ref = ref,
                outbase = outbase
        }
    }

    output {
        Array[File] bam = align.bam
    }
}
```

### Example: `map_by` output

For a DAG with `map_by align "sample" xs` where `xs` has columns `sample, fastq1, fastq2`:

```wdl
workflow main {
    input {
        Array[String] sample
        Array[File] fastq1
        Array[File] fastq2
    }

    Array[Pair[String, Pair[File, File]]] align_grouped = collect_by_key(
        zip(sample, zip(fastq1, fastq2))
    )

    scatter (align_g in align_grouped) {
        String sample_val = align_g.left
        Array[File] fastq1_vals = align_g.right.left
        Array[File] fastq2_vals = align_g.right.right

        call align as align {
            input:
                sample = sample_val,
                fastq1 = fastq1_vals,
                fastq2 = fastq2_vals
        }
    }
}
```

---

## 12. Implementation Order

### Phase 1: Scaffolding + basic task definitions ✓

1. Create `swl/transpile/wdl/` package with `__init__.py`, `cli.py`, `__main__.py`.
2. Implement `transpile_dag_dict()` skeleton: validate, iterate steps, call `_validate_supported`.
3. Implement `_task_to_wdl()`: emit a WDL `task` block with `input {}`, `command <<< >>>`, `output {}`, `requirements {}`.
4. Implement type mapping (`_wdl_type`), interpolation (`_interp_to_wdl`).
5. Test: round-trip a simple single-task DAG (like `partial.swl`).

### Phase 2: Workflow block + call wiring ✓

6. Implement `_dag_to_wdl()`: emit the top-level `workflow { }` block.
7. Implement `_binding_to_wdl_expr()`: handle `Input`, `Literal`, `Field` bindings.
8. Implement workflow input declaration.
9. Test: transpile `function.swl` (three tasks in a pipeline).

### Phase 3: Mapped steps (scatter via `map`) ✓

10. Implement `_mapped_step_to_wdl()`: emit `scatter` block around `call`.
11. Implement scatter input indexing vs broadcast input detection.
12. Implement length-expression derivation from table source.
13. Test: transpile `batch.swl` (map align over samples).

### Phase 4: `map_by` via `collect_by_key()` ✓

14. Implement `_mapped_by_step_to_wdl()`: emit `collect_by_key()` + scatter.
15. Implement zip chain construction for multi-column tables.
16. Implement key/value extraction inside scatter body.
17. Test: transpile a simple `map_by` workflow.

### Phase 5: Sub-workflows ✓

18. Implement `_subworkflow_to_wdl()`: recursively transpile `definition.dag`.
19. Implement multi-workflow emission (sub-workflows before parent).
20. Test: transpile `import_partial.swl` (workflow step with sub-workflow).

### Phase 6: Record + struct ✓

21. Implement `_collect_structs()`: walk DAG for record shapes.
22. Implement `_emit_struct()`: emit WDL `struct` definitions.
23. Implement `Record` binding → struct literal emission.
24. Add merge binding rejection with clear error.
25. Implement `_infer_output_type` for all binding kinds.

### Phase 7: Hardening ✓

26. Handle optional types (`String?`, `File?`) with `?` suffix.
27. Handle special characters in task/workflow names.
28. Handle duplicate aliases across steps (append numeric suffix).
29. Add comprehensive validation in `_validate_supported()`.
30. Handle bash `$` vs WDL `~{}` in scripts (command heredoc passes `$` through).
31. Handle empty inputs/outputs.

---

## 13. Validation (`_validate_supported`)

```python
def _validate_supported(dag):
    for step in dag.steps:
        for name, value in step.bindings.items():
            _validate_binding(value, step.id)

    for name, value in dag.outputs.items():
        _validate_output_binding(value, name)


def _validate_binding(value, step_id):
    if isinstance(value, Merge):
        raise ValueError(
            f'WDL transpilation does not support merge bindings '
            f'(step {step_id}). Flatten merges before transpiling.'
        )


def _validate_output_binding(value, name):
    if isinstance(value, Literal):
        raise ValueError(f'WDL does not support literal workflow outputs: {name}')
    if isinstance(value, Merge):
        raise ValueError(f'WDL does not support merge workflow outputs: {name}')
```

---

## 14. Feature Support Summary

| DAG feature | WDL 1.1 | Status |
|---|---|---|---|
| Task step | `task` with `command <<< >>>` | Implemented (Phase 1) |
| Input binding | Call input block argument | Implemented (Phase 2) |
| Step output binding | `call_alias.output_name` | Implemented (Phase 2) |
| Literal binding | Inline literal value | Implemented (Phase 2) |
| Field projection | `.field_name` member access | Implemented (Phase 2) |
| Record merge | Must be flattened before transpilation | Rejected (Phase 6) |
| Record binding | `struct` literal | Implemented (Phase 6) |
| `map` (scatter) | `scatter` block | Implemented (Phase 3) |
| `map_by` | `collect_by_key()` + scatter | Implemented (Phase 4) |
| Sub-workflow | `workflow` block + `call` | Implemented (Phase 5) |
| CPU | `cpu:` in `requirements {}` | Implemented (Phase 1) |
| Memory | `memory:` in `requirements {}` | Implemented (Phase 1) |
| Container | `container:` in `requirements {}` | Implemented (Phase 1) |
| Time | `time_minutes:` in `requirements {}` | Implemented (Phase 7) |
| Output file path | `File name = "~{var}.ext"` | Implemented (Phase 1) |
| Interpolation (var) | `~{var}` in command/output | Implemented (Phase 1) |
| Interpolation (expr) | `~{expr}` passthrough | Implemented (Phase 1) |
| Table binding | Array inputs + scatter indexing | Implemented (Phase 3) |
| Optional types | `Type?` suffix | Planned Phase 7 |
| Workflow doc/desc | `meta {}` / `parameter_meta {}` | Optional |

---

## 15. Key Risks and Mitigations

| Risk | Mitigation |
|---|---|
| **`collect_by_key()` engine support**: Not all WDL engines fully implement `collect_by_key()` (it was added in 1.1). Cromwell and MiniWDL support it; older engines may not. | Document the WDL 1.1 requirement. Provide a fallback `_validate_supported()` flag to reject `map_by` if engine compatibility is a concern. |
| **Nested `Pair` depth**: Tables with many columns produce deeply nested `Pair[..., Pair[..., ...]]` types that are cumbersome. | Cap at 4 levels of nesting; beyond that, use a synthetic grouping task (bash with `jq`). Or use intermediate `Array` flattening. |
| **Scatter inside scatter**: If a mapped step calls a sub-workflow that itself contains scatters, nested scatters may behave differently across engines. | Test on Cromwell and MiniWDL. Avoid deep nesting; flatten at the SWL level if needed. |
| **Script interpolation conflicts**: WDL uses `~{}` in heredoc commands, but SWL bash scripts use `${}` for shell variables. | WDL `command <<< >>>` treats `~{}` as WDL expressions and passes `$` through literally. Shell `${}` variables are preserved. |
| **Array index bounds**: `scatter (i in range(n))` where `n` is derived from one column assumes all table columns have equal length. | Rely on SWL compiler's compile-time table length validation (spec §Types). |
| **Broadcast scalar inputs in scatters**: Scalars referenced inside a scatter are broadcast automatically by WDL. | Use scalar references directly; no need to wrap in arrays. |
| **Sub-workflow naming**: Multiple sub-workflows in the same file need unique names. | Append numeric suffix for duplicate names. Use `_wf_name()` with collision detection. |
| **Struct name collisions**: Generated struct names from `_struct_name_for()` may collide across structurally similar but semantically distinct records. | Append a hash suffix to the struct name, or use a counter-based naming scheme. |
| **`command <<< >>>` vs `command {}`**: The heredoc form avoids confusion with bash `$`. | Always use `command <<< >>>` with `~{}` placeholders. |
