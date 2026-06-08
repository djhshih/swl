# CWL Transpiler: Implementation Status

## Overview

The CWL transpiler (`python/swl/transpile/cwl/`) converts compiled SWL DAGs into packed CWL v1.0 documents (`cwlVersion: v1.0`) with a `$graph` containing one `Workflow` node and one `CommandLineTool`/`ExpressionTool`/`Workflow` per distinct step. It follows the same architecture as the Snakemake, Nextflow, and WDL transpilers. The transpiler is ~641 lines of Python (`emit.py`).

**Status:** Implemented. Covers simple tasks, pipelines, named ports, resources (cpu/memory/time/image), output path templates, scattered mapped steps, `map_by` (grouped scatter) via ExpressionTool grouping + wrapper CommandLineTool, record bindings via ExpressionTool, nested field projections via `valueFrom`, `expr` interpolation via `InlineJavascriptRequirement`, optional types (`['null', Type]`), sub-workflow inlining, array types (`[file]`, `[str]`, `[int]`, `[float]`), and shell variable interpolation via `$(...)` CWL expressions.

---

## Architecture

```
transpile/cwl/
├── __init__.py     # exports transpile_dag_dict, transpile_dag_file
├── __main__.py     # python -m entry point
├── cli.py          # CLI wrapper (uses transpile/_cli.py)
└── emit.py         # transpiler
```

### Entry points

- `transpile_dag_file(path)` — read DAG JSON from disk, transpile to packed CWL
- `transpile_dag_dict(data, workflow_id='main')` — main transpilation from parsed DAG data

### Design

- **Variable reference forms**: Both `$var` (unbraced) and `${var}` (braced) are handled identically — they are equivalent bash syntax for simple variable names. The `interp_script()` regex matches both forms. `${expr}` (non-word content inside braces) is treated as an expression.
- **Scope resolution**: Each variable reference is classified via `classify_var()`: known inputs → `inputs.name` (`.path` for files), run params → inlined literals, shell builtins (`$HOME`, `$PATH`) → verbatim pass-through, unknown → verbatim pass-through.
- **Shell built-in allowlist**: `_BUILTIN_VARS` in `bash.py` prevents shell builtins from being incorrectly resolved as SWL variables.
- **`InlineJavascriptRequirement`**: Automatically added when the body contains any `$` reference, since the `entry` uses `$(...)` JS expression syntax.

### Key internal functions

| Function | Purpose |
|----------|---------|
| `_tool_to_cwl(step)` | Convert a task step to `CommandLineTool` with `InitialWorkDirRequirement`, `ResourceRequirement`, `DockerRequirement`, `TimeLimit` hints. Recursively packs sub-workflow definitions |
| `_step_to_cwl(workflow_id, step, tool_id, record_map)` | Emit a workflow step entry with `in` array wiring, `scatter` + `scatterMethod: dotproduct` for mapped steps |
| `_workflow_input_to_cwl(workflow_id, name, spec)` | Emit `Workflow.inputs[]` entries with type and optional doc |
| `_workflow_output_to_cwl(workflow_id, name, output)` | Emit `Workflow.outputs[]` with `outputSource` |
| `_step_input_to_cwl(task_id, name, value)` | Convert a DAG binding to a step `in` entry with source and optional `valueFrom` for field projections |
| `_tool_input_to_cwl(tool_id, name, spec)` | Emit tool input port with type mapping |
| `_tool_output_to_cwl(tool_id, name, spec)` | Emit tool output with `outputBinding.glob` as `$(...)` CWL expression |
| `_resource_requirement(run)` | Map `run.cpu` → `coresMin`, `run.memory` → `ramMin` (MB) |
| `_docker_requirement(run)` | Map `run.image` → `dockerPull` |
| `_hints_from_run(run)` | Map `run.time` → `TimeLimit` hint |
| `_binding_source(value)` | Resolve a DAG binding to a CWL `source` string (Input → `#main/name`, StepCall → `#main/step_id`, Field → `#main/root/name`) |
| `_interp_to_cwl_glob(value, input_names=())` | Render a SWL `word` interpolation dict to a CWL expression `$(inputs.var + '.ext')`. Resolves inner variable names in `expr` parts using `input_names` |
| `_has_expr_interpolation(spec)` | Check if an output spec contains `expr` interpolation parts; if so, `InlineJavascriptRequirement` is added |
| `_interpolate_shell(body, step)` | Resolve bash `${var}` / `${expr}` references to a `$(...)` CWL JS expression. Input vars → `inputs.name` / `inputs.name.path`, run params → inlined literals, builtins → verbatim |
| `_validate_supported(dag)` | Validate map_by preconditions, reject literal workflow outputs, validate output interpolation |
| `_emit_record_tool(step_id, binding_name, record, dag)` | Generate an `ExpressionTool` + step entry for a Record binding — constructs JSON from field inputs |
| `_emit_map_by_graph(step, dag)` | Generate ExpressionTool (grouping) + CommandLineTool (wrapper) + two step entries for `map_by` grouped scatter |
| `_gen_wrapper_script(col_names, out_names, original_body)` | Generate the Python wrapper bash script for map_by processing — reads group JSON, iterates rows, calls inner script per row |

---

## SWL → CWL v1.0 Mapping

### Well Supported

| SWL Concept | CWL v1.0 Output | Details |
|---|---|---|
| **Task step** | `CommandLineTool` with `InitialWorkDirRequirement` | Script body materialized as `script.sh` via `$(...)` JS expression |
| **Shell variable interpolation (`${var}`)** | `$(...)` CWL expression in `entry` | Input vars → `inputs.name` / `inputs.name.path`, run params → inlined literals, shell builtins → verbatim. `InlineJavascriptRequirement` added automatically |
| **Shell expression interpolation (`${expr}`)** | `$(...)` with resolved inner vars | `${memory / cpu}` → `(8192 / 2)`. Inner variable names resolved per scope |
| **Named input ports** | Tool `inputs[]` with `type` | `file` → `File`, `str` → `string`, `int` → `int`, `float` → `float` |
| **Named output ports** | Tool `outputs[]` with `outputBinding.glob` | Output paths use `$(inputs.var + '.ext')` CWL expression |
| **`run.cpu`** | `ResourceRequirement.coresMin` | |
| **`run.memory`** | `ResourceRequirement.ramMin` | Value in MB |
| **`run.time`** | `TimeLimit` hint | Emitted as hints entry (not standard ResourceRequirement) |
| **`run.image`** | `DockerRequirement.dockerPull` | |
| **Workflow inputs** | `Workflow.inputs[]` with type and optional doc | |
| **Input binding (workflow input)** | `source: #main/input_name` | Direct workflow input reference |
| **Input binding (literal)** | `default: <value>` | Inline literal value on step input port |
| **Step output binding** | `source: #main/step_id/output_name` | Cross-step reference |
| **Field projection** | `source: #main/root/name` + `valueFrom: $(field_path)` | CWL `valueFrom` with `$(...)` expression for nested access |
| **`word` interpolation in output defaults** | `$(inputs.var + '.literal')` | `word_interp()` from `common.py` with `var_fn=lambda n: f"inputs.{n}"` |
| **`expr` interpolation in output defaults** | `$((inputs.var / 2) + '.literal')` with `InlineJavascriptRequirement` | Expression text resolved with `inputs.` prefix for known input variables; `InlineJavascriptRequirement` added to tool |
| **Workflow output** | `outputs[]` with `outputSource` | Type inferred from source binding |
| **Sub-workflow (no map)** | Embedded `Workflow` node in `$graph` | Recursively transpiled; ID prefixed with step name |
| **Mapped step (scatter)** | `scatter: [ports]` + `scatterMethod: dotproduct` | Workflow inputs for scatter ports promoted to `Array` types |
| **Mapped workflow** | Scattered step referencing sub-workflow | Sub-workflow `Workflow` node in `$graph` + `scatter` step |
| **Mapped lambda** | Generated `Workflow` node + scattered call | Lambda body wrapped in generated sub-workflow, called with scatter |
| **`map_by` (grouped scatter)** | ExpressionTool (grouping) + CommandLineTool (wrapper) | Two-step: groups rows by key via JS, then processes each group via Python wrapper script |
| **Record binding** | `ExpressionTool` with `InlineJavascriptRequirement` | Constructs JSON record from individual field inputs; step entry inserted before consumer |
| **`[file]` / `[str]` / `[int]` / `[float]` array types** | `{type: array, items: File/string/int/float}` | Properly mapped via `to_cwl_type()` |
| **Optional types (`type?`)** | `['null', Type]` union type | Handled via `to_cwl_type()` optional suffix handling |

### Weakly Supported

| SWL Concept | Status | Issue |
|---|---|---|

### Not Supported (explicitly rejected)

| SWL Concept | Behavior |
|---|---|
| **Merge bindings** | Raises `ValueError`: "Unsupported binding for CWL transpilation" (via `_binding_source`) |
| **Literal workflow outputs** | Raises `ValueError`: "literal outputs are not supported" |

### CWL-specific details

- **Packed document format**: All tools and the workflow are flat in `$graph[]`. Tool IDs are `#step_id`, workflow is `#main`. Input/output IDs use hierarchical paths like `#align/fastq1`.
- **Scatter implementation**: All scatter ports are listed explicitly in `scatter: [...]`. `scatterMethod: dotproduct` pairs up array elements by position.
- **`map_by` implementation**: A two-phase emission:
  1. An `ExpressionTool` with `InlineJavascriptRequirement` partitions array inputs by a key column and writes per-group JSON files using `fs.writeFileSync`
  2. A `CommandLineTool` wrapper iterates over group JSON files, calling the original script per row via Python `subprocess.run`
- **Record binding implementation**: An `ExpressionTool` takes record fields as separate inputs, constructs a JSON object via JS, writes it to `record.json`, and returns a `File`. The downstream step references this file.
- **Nested field projections**: `Field(Field(...), ...)` chains are resolved by `field_chain_parts()` to find the root source. The first access is the `source`, remaining accesses become `valueFrom: $(rest.path)`.
- **Optional types**: SWL `type?` suffix is mapped to CWL `['null', Type]` union. This is handled in `to_cwl_type()` for both `File` and scalar types.
- **Requirements emission**: Workflow always declares `ScatterFeatureRequirement` and `MultipleInputFeatureRequirement`. Tools declare `InitialWorkDirRequirement`, `ResourceRequirement`, `DockerRequirement`, `InlineJavascriptRequirement`, and `TimeLimit` hints as needed.

---

## Support Matrix

| Feature | Status | Tests |
|---------|--------|-------|
| Simple task → CommandLineTool | ✅ Well supported | function, pipe, explicit |
| Pipeline chain (|) | ✅ Well supported | function, pipe (equivalence) |
| Named I/O ports | ✅ Well supported | All |
| Resources (cpu/memory/image) | ✅ Well supported | function |
| Time resource hint | ✅ Well supported | function (TimeLimit check) |
| Workflow input/output | ✅ Well supported | All |
| Output path templates (`$(inputs.var)`) | ✅ Well supported | function, pipe, explicit |
| Literal input bindings | ✅ Well supported | partial |
| Field projections | ✅ Well supported | All (source wiring) |
| Sub-workflow (no map) | ✅ Well supported | import_partial |
| Sub-workflow (mapped) | ✅ Supported | panel (scattered sub-workflow) |
| `[file]` / `[str]` / `[int]` / `[float]` types | ✅ Supported | panel (array inputs) |
| Optional types (`type?`) | ✅ Supported | optional_workflow_output test |
| Nested field projections (`a.b.c`) | ✅ Supported | nested_field_binding test |
| Record bindings | ✅ Supported | record_binding_step_input test |
| `expr` interpolation in outputs | ✅ Supported | expr_interpolation_passthrough test |
| Mapped step (scatter) | ✅ Well supported | batch (task, workflow, lambda) |
| `map_by` (grouped scatter) | ✅ Well supported | map_by_transpile_emits, map_by_from_compiled_dag |
| Shell variable interpolation (`${var}`) | ✅ Well supported | shell_interpolation_* tests |
| Shell expression interpolation (`${memory / cpu}`) | ✅ Well supported | shell_interpolation_expression |
| Shell built-in pass-through (`$HOME`) | ✅ Well supported | shell_interpolation_passes_through_builtins |
| Body without `$` (no-op passthrough) | ✅ Supported | shell_interpolation_skips_body_without_dollar, no_dollar_does_not_add_inlinejs |
| Merge bindings | ❌ Rejected | rejects_merged_task_input_binding |
| Literal workflow outputs | ❌ Rejected | (via _validate_supported) |

### Coverage

The test suite (`bash test.sh`):
- Generates `tests/cwl/*.cwl` from all compiled DAGs (`function`, `pipe`, `explicit`, `panel`, `map`)
- Compares `pipe.cwl` and `explicit.cwl` against `function.cwl` for equivalence of golden output
- `panel.cwl` and `map.cwl` are generated but have no golden comparison
- `map_by.cwl` generation is commented out (TODO) but `map_by` is tested via unit tests
- Unit tests: 26 test cases in `tests/unit/swl/transpile/cwl/test_emit.py`

---

## Known Limitations

1. **Merge bindings not supported**: Merge bindings (from `//` operators at the DAG level) are rejected. These should be flattened at compile time before DAG generation.

3. **Literal workflow outputs not supported**: Workflow-level outputs that are literal values (not step outputs) are explicitly rejected.

4. **`map_by` wrapper Python dependency**: The generated wrapper CommandLineTool for `map_by` uses Python 3 with `subprocess` and `json` modules. This introduces a runtime dependency on `python3` that may not be available in all CWL execution environments (e.g., restricted containers).

---

## Example Output

### Simple pipeline (`function.swl`)

Input DAG: `align | sort | call` (3-task pipeline)

Output (`function.cwl`):

```json
{
  "cwlVersion": "v1.0",
  "$graph": [
    {
      "id": "#align",
      "class": "CommandLineTool",
      "baseCommand": ["bash", "script.sh"],
      "inputs": [
        {"id": "#align/fastq1", "type": "File"},
        {"id": "#align/fastq2", "type": "File"},
        {"id": "#align/outbase", "type": "string"},
        {"id": "#align/ref", "type": "File"},
        {"id": "#align/ref_amb", "type": "File"},
        ...
      ],
      "outputs": [{
        "id": "#align/bam",
        "type": "File",
        "outputBinding": {"glob": "$(inputs.outbase + '.bam')"}
      }],
      "requirements": [
        {"class": "InitialWorkDirRequirement",
         "listing": [{"entryname": "script.sh",
                      "entry": "$(\"bwa mem -t \" + 2 + \" \" + inputs.ref.path + \" \" + inputs.fastq1.path + \" \" + inputs.fastq2.path + \" | samtools view -b - > \" + inputs.outbase + \".bam\")"}]},
        {"class": "InlineJavascriptRequirement"},
        {"class": "ResourceRequirement", "coresMin": 2},
        {"class": "DockerRequirement", "dockerPull": "djhshih/seqkit:0.1"}
      ],
      "hints": [{"class": "TimeLimit", "timeLimit": 3150}]
    },
    ...
    {
      "id": "#main",
      "class": "Workflow",
      "inputs": [
        {"id": "#main/fastq1", "type": "File"},
        {"id": "#main/outbase", "type": "string"},
        ...
      ],
      "outputs": [
        {"id": "#main/bam", "type": "File", "outputSource": "#main/sort/bam"},
        {"id": "#main/bai", "type": "File", "outputSource": "#main/sort/bai"},
        {"id": "#main/bcf", "type": "File", "outputSource": "#main/call/bcf"}
      ],
      "requirements": [
        {"class": "ScatterFeatureRequirement"},
        {"class": "MultipleInputFeatureRequirement"}
      ],
      "steps": [
        {"id": "#main/align", "run": "#align",
         "in": [
           {"id": "#main/align/fastq1", "source": "#main/fastq1"},
           {"id": "#main/align/ref", "source": "#main/ref"},
           ...
         ],
         "out": ["#main/align/bam"]},
        {"id": "#main/sort", "run": "#sort", ...},
        {"id": "#main/call", "run": "#call", ...}
      ]
    }
  ]
}
```
