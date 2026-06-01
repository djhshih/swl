# Plan: Decompose `map_by` into `map` at the DAG Level

## Goal

Eliminate `map_by` (and its `MappedStep.group_by` field) as a first-class concept in the DAG JSON. Instead, the compiler decomposes `map_by f key xs` into:

```
_group = group_by(key, xs)    // synthetic grouping step → File[]
f'      = wrap(f)             // synthetic wrapper function: File → rec
result  = map f' _group       // regular map over group JSONs
```

Every transpiler/executor only needs `map` (scatter). The grouping complexity lives in generated bash scripts embedded in the DAG.

---

## Current state

The compiler lowerer produces `ir.Map(function, arg, key=key)` when `key` is set. The forcer creates `MappedStep` with `map.group_by` set:

```python
# force.py:392-397
def _mapped_step_bindings(self, fn, source, key=None):
    info = {'source': _binding_to_public_dict(logical_source)}
    if key is not None:
        info['group_by'] = key
    return {}, info
```

The CWL transpiler rejects any step with `group_by`:

```python
# emit.py:195-196
if getattr(step, 'map', None) is not None and step.map.get('group_by') is not None:
    raise ValueError(...)
```

---

## Design

### Step 1: Generate the grouping step

When the forcer encounters `ir.Map(function, arg, key=key)` (or `MappedStep` building with `key is not None`), instead of setting `map.group_by`, it:

1. Generates a synthetic `StepCall` named `_group_by_{key}_{source}`.
2. The step's script partitions rows by the key column and writes one JSON file per unique key value.
3. The step's inputs are the table columns from `map.source.columns`. Each column becomes an input parameter with its array type (`[file]`, `[str]`, etc.).
4. The step's single output is `_groups: [file]`.

**Generated script structure:**

```bash
#!/bin/bash
# group_by: key=sample, columns=sample,fastq1,fastq2
GROUP_DIR="_groups_$$"
mkdir -p "$GROUP_DIR"

# Read column values from input files or JSON
# ... (see Step 1a below for detailed script)

# Write one JSON per unique key value
for key in $(sorted_unique_keys); do
  cat > "$GROUP_DIR/$key.json" << EOF
{
  "key": "$key",
  "columns": {
    "fastq1": [$(echo "${fastq1_rows[$key]}" | paste -sd,)]
  }
}
EOF
done

echo "$GROUP_DIR"/*.json
```

#### Step 1a: How the grouping script reads its inputs

The table columns reach the DAG as bindings. Each column is either:
- An `Input` (workflow-level table input) — the column values come from the workflow input as an array.
- A `Field` of a `StepCall` (upstream step output) — the column values come from a previous step's array output.

The grouping script needs the actual array values at runtime. Since the columns are arrays, the bindings are wired as array-valued inputs. The generated script reads them from temporary files:

```bash
# For each column input, save to a temp file (one value per line)
# This is necessary because CWL scatter and bash have different array representations
cat "$fastq1_input" > /tmp/fastq1_vals
cat "$fastq2_input" > /tmp/fastq2_vals
cat "$sample_input" > /tmp/sample_vals

# Read into associative arrays
declare -A groups
paste /tmp/sample_vals /tmp/fastq1_vals /tmp/fastq2_vals | while IFS=$'\t' read key fq1 fq2; do
    groups["$key"]+="$fq1,$fq2;"  # append row to group
done
```

**Alternative: single JSON input.** Instead of N separate array inputs, the grouping step could take a single JSON file containing the entire table:

```json
{
  "columns": {
    "sample": ["A", "A", "B"],
    "fastq1": ["A1.fq", "A2.fq", "B1.fq"],
    "fastq2": ["A3.fq", "A4.fq", "B2.fq"]
  },
  "key": "sample"
}
```

This is simpler (one input port instead of N) but requires the predecessor step or workflow input to serialize the table to JSON. The current DAG doesn't have a "serialize table" step, so adding one would mean generating TWO synthetic steps (serialize + group) instead of one. **Rejected — too much indirection.**

**Decision:** N array inputs, one per column. The generated bash reads each from a temp file.

### Step 2: Generate the wrapper function

The wrapper is a synthetic `ir.Function` (kind `task` or `workflow` matching the original). Its body is a generated `ir.Lambda` that:

1. Accepts a single `_group_file: file` input.
2. Reads the group JSON.
3. For each row in the group, extracts per-row scalar values and builds a record for the original function.
4. Calls the original function on each row.
5. Collects outputs.

**Generated script (for a task wrapper):**

```bash
#!/bin/bash
# wrap(original_task): reads group JSON, calls original script per row
GROUP_FILE="$1"

# Read the group
KEY=$(jq -r '.key' "$GROUP_FILE")
ROWS=$(jq -c '.rows[]' "$GROUP_FILE")

OUTPUTS=()
echo "$ROWS" | while read -r row; do
    # Extract per-row values
    FASTQ1=$(echo "$row" | jq -r '.fastq1')
    FASTQ2=$(echo "$row" | jq -r '.fastq2')
    
    # Create temp dir per row
    ROW_DIR="_row_$$"
    mkdir -p "$ROW_DIR"
    
    # Call original tool script
    # (the original script is embedded)
    cd "$ROW_DIR"
    bash -c "$ORIGINAL_SCRIPT" --fastq1 "$FASTQ1" --fastq2 "$FASTQ2"
    
    # Collect output
    # (output names and paths are known from the original task definition)
    cp bam "$ROW_DIR/../${KEY}_bam"
    OUTPUTS+=("$ROW_DIR/../${KEY}_bam")
done

# Output: each collected output file
echo "${OUTPUTS[@]}"
```

**Alternative: emit per-row outputs.** The wrapper could emit each row's output as a separate file for the CWL scatter to collect. But since the wrapper is run once per group (via `map`), the group's results are already collected by the wrapper.

**Decision:** The wrapper produces one output file per output field, aggregated from all rows in the group. Since the wrapper is called once per unique key value (via `map` over group JSONs), the `map` step collects these into per-output arrays.

### Step 3: Replace `MappedStep.group_by` with two regular steps

In `force.py:_emit_mapped_step`, when `key is not None`:

```python
def _emit_mapped_step(self, fn, source, key=None):
    if key is not None:
        return self._emit_map_by(fn, source, key)
    # ... existing map logic ...
```

The new `_emit_map_by` method:

```python
def _emit_map_by(self, fn, source, key):
    # 1. Generate grouping step
    group_step = self._emit_grouping_step(fn, source, key)
    
    # 2. Generate wrapper function
    wrapper_fn = self._build_wrapper_function(fn)
    
    # 3. Emit regular MappedStep over the grouping step's output
    group_source = Field(group_step, '_groups')
    return self._emit_mapped_step(wrapper_fn, group_source, key=None)
```

The wrapper function (`_build_wrapper_function`) creates a synthetic `ir.Function` with:
- Generated name: `_wrap_{original_name}`
- Kind: same as original (`task` or `workflow`)
- Signature: single input `_group_file: file`, outputs = original outputs
- Path: same as original (the wrapper definition is synthetic, but `path` is provenance)
- Body: a generated `ir.Lambda` that reads the group JSON and calls the original

For the wrapper's `body`, when the original is a **task**, the wrapper body is a synthetic `ir.Lambda` whose block contains an `ir.Apply` of the original `ir.Function` to a per-row record, wrapped in a loop. Since the IR is tree-shaped, not imperative, the loop is expressed as a generated bash script in the wrapper's `task` definition rather than in the IR itself.

The wrapper is therefore emitted as a `StepCall` whose `task` field contains a generated `CommandLineTool` definition:

```python
def _build_wrapper_tool_definition(self, original_path, original_signature):
    original_task = self._tool_definition(original_path)
    return {
        'doc': f'Wrapped {original_path} for group processing',
        'body': self._generate_wrapper_script(original_task, original_signature),
        'inputs': {'_group_file': {'type': 'file', 'desc': 'JSON file with group data'}},
        'outputs': original_task['outputs'],  # same outputs as original
        'run': original_task.get('run', {}),
    }
```

### Step 4: Generated grouping script (detailed)

```python
def _generate_group_script(self, columns, key_name):
    """Generate bash script that groups rows by key and writes per-group JSONs."""
    col_names = list(columns.keys())
    col_refs = ' '.join(f'${{{n}_input}}' for n in col_names)
    lines = [
        '#!/bin/bash',
        f'# group_by: key={key_name}, columns={",".join(col_names)}',
        'set -euo pipefail',
        '',
        f'GROUP_DIR="_groups_$$"',
        'mkdir -p "$GROUP_DIR"',
        '',
        '# Read column values into temp files',
    ]
    for name in col_names:
        lines.append(f'{name}_tmp=$(mktemp)')
        lines.append(f'cat "${{{name}_input}}" > "${{{name}_tmp}}"')
    
    lines += [
        '',
        '# Determine unique keys',
        f'cut -f1 "${{{key_name}_tmp}}" | sort -u > "$GROUP_DIR/_unique_keys"',
        '',
        '# Build per-key column value lists',
        f'paste {" ".join(f"${{{n}_tmp}}" for n in col_names)} | while IFS=$\'\\t\' read -r {" ".join(col_names)}; do',
        f'    echo "${{{key_name}}}" >> "$GROUP_DIR/${{{key_name}}}_values"',
        '    for col in ' + ' '.join(col_names) + '; do',
        '        echo "${!col}" >> "$GROUP_DIR/${key}_${col}_values"',
        '    done',
        'done',
        '',
        '# Write per-group JSON files',
        'while read -r key; do',
        '    {',
        f'        echo "{{\\"key\\": \\"$key\\", \\"columns\\": {{"',
        '        first=true',
        '        for col in ' + ' '.join(col_names) + '; do',
        '            $first || echo ","',
        '            first=false',
        f'            echo \\"$col\\": [$(paste -d, -s < "$GROUP_DIR/${{key}}_${{col}}_values")]',
        '        done',
        '        echo "}}"',
        '    } > "$GROUP_DIR/$key.json"',
        'done < "$GROUP_DIR/_unique_keys"',
        '',
        '# Output: list of group JSON files',
        'ls "$GROUP_DIR"/*.json',
    ]
    return '\n'.join(lines)
```

### Step 5: Generated wrapper script (detailed)

```python
def _generate_wrapper_script(self, original_task, original_signature):
    """Generate bash script that reads a group JSON and calls the original tool per row."""
    input_names = list(original_signature.inputs.keys())
    output_names = list(original_signature.outputs.keys())
    original_body = original_task.get('body', '')
    
    lines = [
        '#!/bin/bash',
        f'# wrap: reads group JSON, calls original per row',
        'set -euo pipefail',
        '',
        'GROUP_FILE="$1"',
        '',
        '# Parse group JSON',
        'ROWS=$(python3 -c "',
        'import json, sys',
        'with open(sys.argv[1]) as f:',
        '    group = json.load(f)',
        'columns = group[\"columns\"]',
        'key = group[\"key\"]',
        'num_rows = len(next(iter(columns.values())))',
        'for i in range(num_rows):',
        '    row = {col: values[i] for col, values in columns.items()}',
        '    print(json.dumps(row))',
        f'" "$GROUP_FILE")',
        '',
    ]
    
    # Collect per-row outputs
    for out_name in output_names:
        lines.append(f'{out_name}_files=()')
    
    lines += [
        '',
        '# Process each row',
        'ROW_INDEX=0',
        'while read -r row_json; do',
        '    ROW_DIR="_row_$$_$ROW_INDEX"',
        '    mkdir -p "$ROW_DIR"',
        '    cd "$ROW_DIR"',
        '',
    ]
    
    # Extract values per input
    for name in input_names:
        lines.append(f'    {name}=$(echo "$row_json" | python3 -c "import json,sys; print(json.load(sys.stdin)[\'{name}\'])")')
    
    lines += [
        '',
        f'    # Call original tool',
        f'    bash -c {repr(original_body)} \\',
    ]
    for name in input_names:
        lines.append(f'        --{name} "${{{name}}}" \\')
    
    lines += [
        '    cd ..',
    ]
    
    for out_name in output_names:
        lines.append(f'    {out_name}_files+=("$ROW_DIR/{out_name}")')
    
    lines += [
        '    ROW_INDEX=$((ROW_INDEX + 1))',
        'done <<< "$ROWS"',
        '',
    ]
    
    # Output concatenation
    lines.append('# Output collected files')
    for out_name in output_names:
        lines.append(f'cat "${{{out_name}_files[@]}}" > "collected_{out_name}"')
        lines.append(f'echo "collected_{out_name}"')
    
    return '\n'.join(lines)
```

---

## Changes to the DAG JSON schema

**Removed:** `map.group_by` field from `MappedStep`. The `group_by` key is no longer a serialized field.

**Unchanged:** Everything else. `MappedStep` still has `map.source`, `input_schema`, `output_schema`. The only difference is that `map.group_by` is never set — the grouping is handled by a preceding `StepCall`.

A `map_by` workflow now produces a DAG with two extra steps:

```json
{
  "version": "1.0",
  "inputs": {
    "sample": {"type": "[str]"},
    "fastq1": {"type": "[file]"},
    "fastq2": {"type": "[file]"}
  },
  "steps": [
    {
      "id": "_group_by_sample",
      "type": "task",
      "script": "#!/bin/bash\n... group rows by sample ...",
      "deps": [],
      "inputs": {
        "sample": {"type": "[str]"},
        "fastq1": {"type": "[file]"},
        "fastq2": {"type": "[file]"}
      },
      "outputs": {
        "_groups": {"type": "[file]"}
      },
      "bindings": {
        "sample": {"source": "input", "name": "sample"},
        "fastq1": {"source": "input", "name": "fastq1"},
        "fastq2": {"source": "input", "name": "fastq2"}
      }
    },
    {
      "id": "_wrap_align",
      "type": "task",
      "script": "#!/bin/bash\n... read group JSON, call align per row ...",
      "deps": ["_group_by_sample"],
      "inputs": {
        "_group_file": {"type": "file"}
      },
      "outputs": {
        "_wrap_out": {"type": "file"}
      },
      "bindings": {
        "_group_file": {"step": "_group_by_sample", "output": "_groups"}
      },
      "map": {
        "source": {"source": "step", "step": "_group_by_sample", "output": "_groups"}
      }
    }
  ],
  "outputs": {
    "_wrap_out": {"step": "_wrap_align", "output": "_wrap_out"}
  }
}
```

---

## Changes to the forcer (`force.py`)

### New method: `_emit_map_by`

```python
def _emit_map_by(self, fn, source, key):
    column_bindings = self._extract_column_bindings(source)
    group_step = self._emit_grouping_step(key, column_bindings, fn)
    wrapper_fn = self._build_wrapper_function(fn)
    group_source = Field(group_step, '_groups')
    return self._emit_mapped_step(wrapper_fn, group_source, key=None)
```

### New method: `_emit_grouping_step`

```python
def _emit_grouping_step(self, key, column_bindings, original_fn):
    col_names = sorted(column_bindings.keys())
    call_id = self._step_id(f'_group_by_{key}')
    inputs = {}
    for name, binding in column_bindings.items():
        col_type = original_fn.signature.inputs.get(name)
        inputs[name] = {
            'type': f'[{col_type.type.value}]' if col_type and col_type.type else '[str]',
            'desc': f'Column {name} for grouping by {key}'
        }
    step = StepCall(
        id=call_id,
        path='_group_by',
        bindings=dict(column_bindings),
        outputs=['_groups'],
        run={},
        task={
            'doc': f'Group rows by {key}',
            'body': self._generate_group_script(column_bindings, key),
            'inputs': inputs,
            'outputs': {'_groups': {'type': '[file]', 'desc': f'Per-group JSON files'}},
            'run': {},
        },
        deps=sorted(self._step_dependencies(column_bindings)),
    )
    self.steps.append(step)
    return step
```

### New method: `_build_wrapper_function`

```python
def _build_wrapper_function(self, fn):
    original = fn.function
    wrapper_name = f'_wrap_{original.name}'
    wrapper_signature = TaskSignature(
        inputs={'group_file': Param(['group_file'], TypeKind.FILE)},
        outputs=dict(original.signature.outputs),
        run={},
    )
    wrapper = ir.Function(
        name=wrapper_name,
        kind='task',
        signature=wrapper_signature,
        path=original.path,
        is_batch=False,
        body=ir.Lambda(
            param='_group_file',
            body=ir.Block([], ir.Name('_group_file')),
        ),
    )
    return ForcedFunction(wrapper, None, wrapper_signature)
```

The wrapper's `body` is a dummy `ir.Lambda` because the actual work happens in the generated bash script stored in the `task` definition. The `ForcedFunction` wrapper is emitted via `_emit_mapped_step` → `_emit_mapped_step` calls `_tool_definition(path)` which would try to re-read the original source file. To avoid this, `_emit_mapped_step` must use the generated tool definition instead of re-reading the file.

**Required change in `_emit_mapped_step`:** When the function is a wrapper (detected by `fn.function.name.startswith('_wrap_')`), use the generated tool definition from `_build_wrapper_tool_definition()` instead of calling `_tool_definition(path)`.

### Modified: `_emit_mapped_step`

```python
def _emit_mapped_step(self, fn, source, key=None):
    target = fn.function
    outputs = list(target.signature.outputs.keys())
    step_id = self._step_id(target.name)
    
    if key is not None:
        return self._emit_map_by(fn, source, key)
    
    if target.name.startswith('_wrap_'):
        task_def = self._build_wrapper_tool_definition(target)
    else:
        task_def = self._tool_definition(target.path) if target.kind == 'task' else self._workflow_definition(target)
    
    # ... rest of existing _emit_mapped_step, but no group_by in map_info ...
    
    return step
```

### Removed from `_mapped_step_bindings`

```python
def _mapped_step_bindings(self, fn, source, key=None):
    logical_source = self._logical_table_source(source)
    info = {'source': _binding_to_public_dict(logical_source)}
    # key is no longer stored in map_info
    return {}, info
```

---

## Changes to the lowerer (`lower.py`)

### Minimal change: `ir.Map` with `key` remains

The lowerer still produces `ir.Map(function, arg, key=key)` for `map_by`. The decomposition happens in the forcer, not the lowerer. This keeps the IR clean (one concept: map, optionally with a key).

The `ir.Map` node's `key` field is preserved as an IR-level concept only. It is consumed by the forcer and never reaches the DAG.

---

## Changes to the CWL transpiler (`emit.py`)

### Removed: `map_by` rejection

```python
# DELETE this block
if getattr(step, 'map', None) is not None and step.map.get('group_by') is not None:
    raise ValueError(...)
```

No `group_by` field exists in the DAG, so this check is never triggered.

### Unchanged: `map` handling

The CWL transpiler's existing `scatter` + `scatterMethod: dotproduct` handling for `MappedStep` works as-is for the new decomposed DAG:

- The grouping step is a regular `StepCall` (no `map` field), so it emits as a normal `CommandLineTool`.
- The mapped wrapper step is a regular `MappedStep` with `map.source = {step: "_group_by_sample", output: "_groups"}`, which is a step_output binding. The transpiler scatters over this `File[]` output.

The CWL transpiler needs no new logic — it already handles scattered steps that iterate over a `File[]` from a previous step's output.

---

## Edge cases and risks

| Edge case | Handling |
|-----------|----------|
| **Empty table** (no rows) | The grouping script produces zero JSON files. `map` over an empty array produces zero outputs. The resulting DAG has no outputs for the wrapper step. The workflow output would be an empty array, which is correct. |
| **Single group** (all rows have same key) | One JSON file written. `map` runs the wrapper once. Correct. |
| **Many groups** (thousands of unique keys) | Thousands of JSON files. CWL scatter handles this, but file system pressure and `jq`/`python3` overhead per group may be a concern. The per-group overhead is bounded by the wrapper script's startup cost. |
| **Key column contains special characters** (`/`, `\`, etc.) | The grouping script must sanitize key values for use as filenames. Use `md5sum` or base64-encode the key for the filename. |
| **Original function has array-typed inputs** (e.g., `merge.sh: bam [file]`) | In `map_by` context, these become arrays-of-arrays. The wrapper must handle this. For initial implementation, reject this case in `_emit_map_by` — detect array-typed inputs and raise a clear error. |
| **Original function has array-typed outputs** | No issue — the wrapper collects per-row outputs and concatenates them. Each group produces one output per output field, which `map` collects into an array. |
| **Workflow function (not task)** | The wrapper is generated as a synthetic `Workflow` step with an embedded `WorkflowDefinition`. The grouping step is unchanged. The wrapper's `definition` field contains a generated sub-workflow that: (1) reads the group JSON via an inline task, (2) calls the original sub-workflow per row. |
| **`expr` interpolation in original function** | The wrapper script inherits the original script body. If the original uses `${cpu * 2}`, the wrapper passes it through. The CWL transpiler already rejects `expr` in output paths; no new rejection logic needed. |
| **Input column types** | The grouping script needs to know column types (file vs str vs int) to write correct JSON. The column types come from the table source's column bindings. Files should be serialized as `{"class": "File", "path": "/path/to/file"}` for CWL compatibility. |
| **No unique sort order** | The grouping script uses `sort -u` which produces lexicographic order. Group order in the output may not match input order. This matches the spec (order not guaranteed). |
| **`map_by` with partial application** (`map_by f "key"` where `f` is already partially bound) | The wrapper must carry the pre-bound inputs into each per-row invocation. The `_build_wrapper_function` needs access to `fn.bound` and must embed the bound values in the generated wrapper script. |

---

## Implementation order

1. **Add `_generate_group_script` and `_generate_wrapper_script`** to `force.py` as standalone functions (testable without the Forcer class).
2. **Add `_emit_map_by`, `_emit_grouping_step`, `_build_wrapper_function`, `_build_wrapper_tool_definition`** to the `Forcer` class.
3. **Modify `_emit_mapped_step`** to call `_emit_map_by` when `key is not None`.
4. **Modify `_mapped_step_bindings`** to not store `group_by`.
5. **Remove `group_by` from `MappedStep`** (or keep the field but never set it — simpler to keep it for forward compatibility).
6. **Update `dag.py` serialization** if `MappedStep.group_by` is removed from the dataclass.
7. **Remove the rejection from `emit.py:_validate_supported`** (the DAG no longer has `group_by` to reject).
8. **Add tests**: `map_by` produces a DAG with `_group_by_*` step + scattered `_wrap_*` step; round-trip serialization; CWL transpilation produces valid CWL with scatter.

## Open questions

1. **Script readability**: The generated bash scripts are complex. Should we use Python scripts instead (via `python3 -c` embedded in bash)? Python is more reliable for JSON handling and array manipulation. **Tentative:** Use `python3` for JSON parsing in the wrapper, plain bash + `jq` for the grouping step (since `jq` is standard in CWL environments).

2. **File serialization format**: How should `file`-typed columns be serialized in the group JSON? CWL's `File` type includes metadata (`class`, `path`, `size`, `checksum`). For the grouping JSON, should we use `{"class": "File", "path": value}` for CWL compatibility, or just the path string? **Decision:** Use path strings in the group JSON (simpler). The wrapper script can add `{"class": "File"}` wrapping if the target requires it.

3. **Binding column extraction**: `_extract_column_bindings(source)` needs to walk the source binding (table → columns) and produce a flat dict of column name → binding. This logic partially exists in `_logical_table_source`. Should be straightforward to extract.

4. **Step deduplication**: If the same table is used by multiple `map_by` calls with different keys (e.g., `map_by f "sample" xs` and `map_by g "region" xs`), the grouping step is duplicated. Should we cache grouping steps by (source key, column signature)? For V1, no dedup — the extra step is small. Add dedup via `_step_call_key` if needed.
