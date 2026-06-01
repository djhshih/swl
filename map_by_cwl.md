# Supporting `map_by` in CWL Transpilation

## What `map_by` means

In SWL, `map_by f key xs` partitions the logical rows of table `xs` by equality of values in column `key`, applies function `f` to each partition (group), and produces one output row per group:

```
xs : tab{ sample: [str], fastq1: [file], fastq2: [file] }
f  : rec{ sample: [str], fastq1: [file], fastq2: [file] } -> rec{ bam: file }

map_by f "sample" xs
  → groups rows by unique sample value
  → applies f to each group (f receives array-typed fields: the entire group's data)
  → output tab{ bam: [file] } — one row per group
```

This is fundamentally different from `map`:
- `map f xs`: `f` receives one **row** at a time (scalar fields). CWL `scatter` naturally expresses this.
- `map_by f key xs`: `f` receives one **group** at a time (array fields). CWL has no native grouped scatter.

## How `map_by` appears in the DAG JSON

The compiler records `map_by` as a `MappedStep` with:

```json
{
  "type": "mapped_task",
  "map": {
    "source": { "source": "table", "name": "table", "columns": {"sample": {"source": "input", "name": "sample"}, ...} },
    "group_by": "sample"
  },
  "input_schema": { "sample": "str", "fastq1": "file", "fastq2": "file" },
  "output_schema": { "bam": "file" },
  ...
}
```

The `input_schema` and `output_schema` contain **per-row types** (the same types the function uses in a non-batched `map`). The `group_by` field signals that rows must be grouped before the function is applied. The `map.source` identifies the table being grouped.

The current CWL transpiler rejects any step with `group_by` set.

## Why CWL can't express `map_by` natively

CWL `scatter` (v1.0–v1.2) supports three `scatterMethod` values:

| `scatterMethod` | Behavior |
|----------------|----------|
| `dotproduct` | Takes N arrays, passes element `i` from each to invocation `i` |
| `nested_crossproduct` | Takes N arrays, invokes every combination |
| `flat_crossproduct` | Like nested but flattens one level |

None of these do grouping. They all assume the input data is already partitioned into per-invocation elements. There is no CWL built-in for "group rows by column value, then scatter over groups."

## Approach: Pre-grouping step + scatter

The only way to implement `map_by` in CWL is to split it into two phases:

1. **Phase 1 — Grouping**: A step that partitions the table rows by the key column and emits per-group data.
2. **Phase 2 — Scatter**: A scattered step that runs the mapped function on each group.

### Phase 1: The grouping step

Use a CWL `ExpressionTool` with `InlineJavascriptRequirement`. It takes the table columns (as arrays) and the key name, groups rows by the key, and writes one JSON file per group.

```
input:  sample: string[],  fastq1: File[],  fastq2: File[]
output: groups: File[]     (one JSON file per unique key value)
```

Each JSON file contains the group's data:

```json
{
  "key": "sample_A",
  "columns": {
    "fastq1": [{"class": "File", "path": "/paths/A1.fq"}, {"class": "File", "path": "/paths/A2.fq"}],
    "fastq2": [{"class": "File", "path": "/paths/A3.fq"}, {"class": "File", "path": "/paths/A4.fq"}]
  }
}
```

The ExpressionTool is generic: it needs to know only the column bindings (available from `map.source.columns` in the DAG) and the key name (`map.group_by`).

### Phase 2: The scattered step

Scatter over the group JSON files. For each group file, the mapped function needs to process the group's data.

This requires a **generated wrapper CommandLineTool** that:
1. Reads a group JSON file
2. Converts the array-typed group columns into the per-row arguments expected by the original tool
3. Calls the original tool once per row (or however the function processes the group)
4. Collects the outputs

The wrapper tool is specific to the mapped function. It is generated at transpile time by combining:
- The original tool's definition (from `step.task` in the DAG)
- The group JSON schema (from `map.source.columns`)
- The output schema (from `step.output_schema`)

## Design option

### ExpressionTool grouping + scattered wrapper

**Flow:**
```
[table columns] → GroupBy (ExpressionTool) → [group JSONs] → ScatteredWrapper → [results]
                                                              (scatter over JSONs)
```

**Pros:**
- Uses standard CWL `scatter` for parallelism
- Grouping logic is isolated in a reusable ExpressionTool
- The wrapper tool is similar to what the compiler already generates for `map` + lambda

**Cons:**
- Requires `InlineJavascriptRequirement` (widely supported but not universal)
- The wrapper tool adds complexity to the CWL output
- Group JSON files need to be written and read (I/O overhead)

**CWL output structure:**
```yaml
class: Workflow
requirements: [ScatterFeatureRequirement, InlineJavascriptRequirement]

inputs:
  sample: string[]
  fastq1: File[]
  fastq2: File[]

steps:
  group_by_sample:
    run:
      class: ExpressionTool
      requirements: [InlineJavascriptRequirement]
      inputs:
        sample: string[]
        fastq1: File[]
        key_name: string
      outputs:
        groups:
          type: { type: array, items: File }
          outputBinding: { glob: "groups/*.json" }
      expression: |
        ${ ... partition rows by key, write per-group JSONs ... }
    in:
      sample: sample
      fastq1: fastq1
      key_name: { default: "sample" }
    out: [groups]

  process_group:
    run: generated_wrapper.cwl
    in: { group_file: groups }
    scatter: [group_file]
    out: [bam]

outputs:
  bam:
    type: File[]
    outputSource: process_group/bam
```

## What needs to change in the transpiler

The current rejection in `emit.py:_validate_supported`:

```python
if getattr(step, 'map', None) is not None and step.map.get('group_by') is not None:
    raise ValueError(...)
```

Would be replaced by logic that:

1. **Extracts metadata**: reads `step.map.group_by`, `step.map.source.columns`, `step.input_schema`, `step.output_schema`, and the function's tool definition from `step.task`.

2. **Generates the grouping ExpressionTool** (Option A):
   - Creates inline CWL with `InlineJavascriptRequirement`
   - The JavaScript expression groups rows by the key column
   - Emits per-group JSON files via `fs.writeFileSync`
   - Outputs a `File[]` of group JSONs with `glob: "groups/*.json"`

3. **Generates the scattered wrapper CommandLineTool**:
   - Takes a single `File` input (the group JSON)
   - Has a generated script that:
     a. Reads the JSON
     b. For each row in the group, calls the original tool's script with the row's values
     c. Collects outputs
   - The wrapper's script can be generated by inspecting `step.task.inputs` and `step.task.outputs` to produce a bash script that loops over rows.

4. **Emits the CWL workflow steps**:
   - Grouping step (ExpressionTool)
   - Scattered processing step (generated wrapper)
   - Scatter on the group file input

5. **Drops unsupported features**: If the group JSON approach can't handle some feature (e.g., `expr` interpolation in the function's arguments), it should still reject with a clear error rather than silently producing incorrect CWL.

## Key risks and edge cases

- **The function uses `expr` interpolation** (e.g., `${outbase / 2}`): The transpiler already rejects `expr` in CWL output paths (`_interp_to_cwl_glob`). This rejection should also apply to `map_by` wrapped steps.

- **The function is a workflow, not a task**: Both need a wrapper to bridge group data to the function's per-row input ports. For a task, the wrapper is a `CommandLineTool` with a generated bash script. For a workflow, the wrapper is a `Workflow` (using `SubworkflowFeatureRequirement`) that reads the group JSON via an inline `ExpressionTool` and calls the inner sub-workflow per row. The transpiler already generates synthetic sub-workflows for `map` + lambda (`lower.py:_generated_callable_from_lambda`). No reason to reject workflow steps.

- **The function's inputs include array-typed parameters** (e.g., `merge.sh` with `bam [file]`): In `map_by` context, these become arrays-of-arrays. The wrapper must handle this. For initial support, reject `map_by` when the function's input schema includes array types (since the column values would be arrays of arrays, which is not representable in the group JSON).

- **Large number of groups**: Each group produces a separate JSON file. With thousands of groups, this creates many small files. The wrapped step should handle this gracefully (CWL scatter can handle thousands of elements, but file system pressure may be a concern).

- **The grouping key column is not in the input schema**: The DAG's `map.source.columns` may include columns not listed in `input_schema`. The ExpressionTool must use `input_schema` to determine which columns to include in the per-group data.
