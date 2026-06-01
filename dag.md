# SWL Compiled DAG JSON Specification

The DAG JSON is the **sole** output of the SWL compiler. It is a self-contained, minimal representation of a workflow as a directed acyclic graph of step invocations. Transpilers and execution engines must work entirely from DAG JSONs — no source files (`.swl`, `.sh`) may be read at transpilation or execution time.

---

## Design Principles

1. **Self-contained**: All task definitions, tool metadata, and scripts are embedded. A DAG JSON can be transpiled or executed without access to the original source files.
2. **Transpiler-neutral**: The schema represents function application, record merging, field projection, and batch semantics in a target-agnostic way.
3. **Minimal but sufficient**: Every field exists to serve at least one transpilation target or executor. Target-specific enrichment (e.g., CWL `scatterMethod`, Nextflow `fromChannel`) is added by the transpiler, never by the schema.

---

## 1. Top-Level Structure

```json
{
  "version": "1.0",
  "inputs": { "<name>": <InputSpec>, ... },
  "steps": [ <StepSpec>, ... ],
  "outputs": { "<name>": <Binding>, ... }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `version` | string | yes | Schema version, currently `"1.0"` |
| `inputs` | object | yes | Map of workflow input name → InputSpec |
| `steps` | array | yes | Ordered list of step specifications |
| `outputs` | object | yes | Map of workflow output name → Binding |

---

## 2. Data Model

The DAG JSON's data model represents all values that flow through the workflow graph.

### 2.1 Scalar Types

The type system uses SWL type names. These are normalized — `"File"` and `"string"` from task annotations are converted to `"file"` and `"str"`.

| Type string | JSON type | Description |
|-------------|-----------|-------------|
| `"file"` | string | File path |
| `"str"` | string | String |
| `"int"` | integer | Integer |
| `"float"` | number | Floating point |
| `"[file]"` | array | Array of file paths (table column) |
| `"[str]"` | array | Array of strings |
| `"[int]"` | array | Array of integers |
| `"[float]"` | array | Array of floats |
| `null` | null | Unknown/inferred type |

When a `?` (optional) annotation is present, the type string is unchanged in the DAG; the optionality is a validational concern that transpilers may inspect via a separate flag if needed. Currently optionality is preserved only at the semantic level, not serialized in the DAG.

### 2.2 InputSpec

```json
{
  "type": "file" | "str" | "int" | "float" | "[file]" | "[str]" | "[int]" | "[float]" | null,
  "desc": "description string or null"
}
```

Declares a named workflow input.

### 2.3 ParamSpec

Describes an input or output parameter of a task or workflow step.

```json
{
  "type": "file" | "str" | "int" | "float" | "[file]" | "[str]" | "[int]" | "[float]" | null,
  "desc": "human-readable description or null",
  "default": <Interpolation> or null
}
```

The `default` field is an Interpolation (see §6) representing the default value expression. For output parameters, this is the expected file path or glob pattern.

### 2.4 RunParamSpec

Runtime resource parameter. Resources are always resolved to concrete values by compile time.

```json
{
  "type": "memory" | "time" | "int" | "str" | null,
  "value": 8192,
  "desc": "description or null"
}
```

| `type` | `value` semantics | Units |
|--------|-------------------|-------|
| `"memory"` | Memory in MB | MB (integer) |
| `"time"` | Wall time in minutes | Minutes (integer) |
| `"int"` | Core count | Cores (integer) |
| `"str"` | Docker image | Pull string (e.g., `"ubuntu:22.04"`) |

### 2.5 Binding

A Binding is a recursive data structure that describes how a value is wired.

**`input`**: Reference to a workflow input.
```json
{ "source": "input", "name": "input_name" }
```

**`step_output`**: Reference to an output of a previous step.
```json
{ "step": "step_id", "output": "output_field_name" }
```

**`literal`**: A constant value.
```json
{ "source": "literal", "value": <any JSON value> }
```

**`field`**: Field projection from another binding expression.
```json
{ "source": "field", "field": "field_name", "value": <Binding> }
```
Used when a step input requires `record.field` rather than the whole record.

**`merge`**: Record merge, right-overrides-left (`r1 // r2`).
```json
{ "source": "merge", "left": <Binding>, "right": <Binding> }
```
Used by pipeline desugaring (`A | B | C`) to accumulate outputs. Transpilers should flatten merge trees to flat record structures where the target language supports it.

**`record`**: A constructed record literal.
```json
{ "source": "record", "fields": { "<name>": <Binding>, ... } }
```
Used when partial application creates a concrete input record or when a lambda produces a result record.

**`table`**: A batch table source — a collection of named column bindings.
```json
{ "source": "table", "name": "table_name", "columns": { "<name>": <Binding>, ... } }
```
Used in mapped steps to describe the table being scattered over.

### 2.6 Feature: Data Model Implementation Status

| Sub-feature | Implemented | Tested |
|-------------|-------------|--------|
| Scalar types (file/str/int/float) in InputSpec and ParamSpec | Yes | `test_transpile_function_workflow` |
| Array types ([file]/[str]/[int]/[float]) | Yes | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` |
| Input type/desc metadata from step specs | Partial: types inferred from step inputs, but inputs without matching step input names lack type/desc | implicit |
| RunParamSpec with resolved values (cpu/memory/time/image) | Yes (cpu, memory, image, time) | `test_transpile_function_workflow` |
| Input binding | Yes | all CWL tests |
| Step output binding | Yes | all CWL tests |
| Literal binding | Yes | `test_rejects_non_task_workflow_output` |
| Field binding | Yes (input → field, step → field) | partial: nested fields not transpilable |
| Merge binding | Yes (serialized) but no target supports it | CWL transpiler rejects merge bindings |
| Record binding | Yes (serialized) but no target supports it | CWL transpiler rejects record bindings |
| Table binding (batch source) | Yes | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` |
| Optional type (`file?`, `str?`) serialization | No | not tested |

---

## 3. Computation Model: Steps

A step is a node in the workflow DAG. There are four step types, each with a distinct `type` discriminator:

| `type` value | Meaning |
|--------------|---------|
| `"task"` | Runs a bash script (in-process execution) |
| `"workflow"` | Invokes a sub-workflow (inlines or calls sub-DAG) |
| `"mapped_task"` | A task step wrapped in `map` or `map_by` (batch over table rows) |
| `"mapped_workflow"` | A workflow step wrapped in `map` or `map_by` |

All step types share a common base structure:

```json
{
  "id": "align",
  "type": "task",
  "script": "bwa mem ...\nsamtools sort ...\n",
  "deps": [],
  "inputs": { "fastq1": { "type": "file", "desc": "read 1" }, ... },
  "outputs": { "bam": { "type": "file", "default": { "kind": "word", "parts": [...] } }, ... },
  "run": { "cpu": { "type": "int", "value": 4 }, ... },
  "bindings": { "fastq1": { "source": "input", "name": "fastq1" }, ... }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | yes | Unique step identifier within the DAG |
| `type` | string | yes | Step type discriminator |
| `script` | string | yes | Full bash script body (task) or `""` (workflow) |
| `deps` | array of strings | yes | Step IDs this step depends on (may be empty). Derived from bindings — any binding that references a step's output adds that step as a dependency. |
| `inputs` | object of ParamSpec | yes | Parameter declarations for tool inputs |
| `outputs` | object of ParamSpec | yes | Parameter declarations for tool outputs |
| `run` | object of RunParamSpec | yes | Runtime resource declarations |
| `bindings` | object of Binding | yes | Wires every named input parameter to a value. Keys match `inputs`. |

### 3.2 Mapped Steps (batch semantics)

A mapped step adds a `map` field to the step object:

```json
{
  "id": "align_all",
  "type": "mapped_task",
  "map": {
    "source": { "source": "input", "name": "samples" },
    "group_by": "sample_id"
  },
  "input_schema": { "fastq1": "file", "fastq2": "file" },
  "output_schema": { "bam": "file" },
  ...
}
```

| Field | Type | Description |
|-------|------|-------------|
| `map.source` | Binding | The table source. Typically an `input` (single table input) or a `table` binding with named columns. |
| `map.group_by` | string or absent | If present, this is a `map_by` — rows are grouped by this column before the function is applied. If absent, this is a simple `map` — the function applies independently to each row. |

The `input_schema` and `output_schema` contain the **per-row** types (not array-wrapped). Transpilers use these to:
- Declare scatter port types (CWL: `type: File` not `type: {type: array, items: File}`)
- Construct per-row input records for the mapped function
- Derive output table column types

### 3.3 Workflow Steps (sub-workflow invocation)

A workflow step adds a `definition` field:

```json
{
  "id": "pipeline",
  "type": "workflow",
  "definition": {
    "class": "Workflow",
    "dag": { "version": "1.0", "inputs": {...}, "steps": [...], "outputs": {...} },
    "inputs": { ... },
    "outputs": { ... },
    "doc": "sub-workflow doc"
  },
  ...
}
```

The `definition` recursively contains:
- `class`: always `"Workflow"`
- `dag`: the full DAG structure (same schema, recursively self-contained)
- `inputs`/`outputs`: the sub-workflow's own interface parameter specs
- `doc`: the sub-workflow's documentation string

The top-level `inputs`, `outputs`, `script`, and `run` fields are still present and duplicate the interface-level info from `definition` for convenience. Transpilers that inline sub-workflows (CWL) use the `definition.dag`; transpilers that call sub-workflows as separate entities (some WDL patterns) use the top-level interface fields.

### 3.4 Feature: Step Model Implementation Status

| Sub-feature | Implemented | Tested |
|-------------|-------------|--------|
| Task step (type: "task") | Yes | `test_transpile_function_workflow` |
| Workflow step with recursive embedded DAG | Yes | `test_imported_workflow_output_transpiles_as_workflow_step` |
| Mapped task (type: "mapped_task", simple `map`) | Yes | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` |
| Mapped workflow (type: "mapped_workflow", simple `map`) | Yes | `test_batch_mapped_workflow_emits_scattered_subworkflow` |
| Generated sub-workflow for map + lambda | Yes | `test_batch_mapped_simple_lambda_emits_generated_scattered_subworkflow` |
| Root-level partial map (partial application of `map`) | Yes | `test_root_partial_map_transpiles_as_scattered_subworkflow` |
| Mapped step with `map_by` (group_by) | Yes (compiler produces group_by) but no transpiler supports it | `test_map_by_transpile_reports_explicit_grouping_not_supported` |
| DAG circularity validation | No | not tested |
| Dependency derivation from bindings | Yes (implicit from step references in bindings) | all CWL tests |
| Partial application with concrete record bindings | Yes | `test_partial_workflow_transpiles` |
| Lambda materialization to synthetic function | Yes | `test_batch_mapped_simple_lambda_emits_generated_scattered_subworkflow` |

---

## 4. Data Flow: Pipeline Semantics

### 4.1 Pipeline (`A | B | C`)

A pipeline desugars at compile time into the equivalent `Apply` + `Update` form (spec §Pipeline), which is then lowered and forced. The DAG never contains a pipeline node — it is always flattened into individual steps with merge bindings:

```
# A | B | C  becomes in the DAG:
#   step A: { outputs: [a1, a2] }
#   step B: { bindings: { ..., a1: {step: "A", output: "a1"}, ... } }
#   step C: { bindings: { ..., a1: {step: "A", output: "a1"}, a2: {step: "B", output: "b1"}, ... } }
#   outputs: { a1: {step: "A", output: "a1"}, a2: {step: "B", output: "b1"}, c1: {step: "C", output: "c1"} }
```

Record merge bindings (`{source: "merge", ...}`) may appear during pipeline desugaring but should be resolved to flat bindings by forcing. If a merge binding does appear in the DAG, it represents a situation where the desugared merge could not be flattened (e.g., merging two non-step sources).

### 4.2 Feature: Pipeline Semantics Implementation Status

| Sub-feature | Implemented | Tested |
|-------------|-------------|--------|
| Chain desugaring (A \| B → lambda block) | Yes | `force.py:_chain_to_lambda_block` |
| Merge flattening in outputs (`a // b // c` → flat record) | Partial: chain desugaring produces flat outputs; hand-written equivalent may not flatten fully | `pipe.swl` vs `function.swl` divergence known |
| Merge flattening in step bindings | No (merge bindings serialized but rejected by CWL transpiler) | `test_rejects_merged_task_input_binding` |

---

## 5. Batch Semantics: `map` and `map_by`

### 5.1 `map f xs`

Desugared to a `MappedStep` with `map.group_by` absent. The `map.source` binding identifies the table input. The `input_schema` and `output_schema` carry per-row types.

### 5.2 `map_by f key xs`

Same structure as `map` but with `map.group_by` set to the grouping key column name. The transpiler must:
1. Group logical rows of the input table by the named column
2. Apply the function to each group
3. Reassemble results

### 5.3 Feature: Batch Semantics Implementation Status

| Sub-feature | Implemented | Tested |
|-------------|-------------|--------|
| `map` on task | Yes | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` |
| `map` on workflow | Yes | `test_batch_mapped_workflow_emits_scattered_subworkflow` |
| `map` on lambda (generated sub-workflow) | Yes | `test_batch_mapped_simple_lambda_emits_generated_scattered_subworkflow` |
| `map` on partial application (root `map f`) | Yes | `test_root_partial_map_transpiles_as_scattered_subworkflow` |
| `map_by` compile-time lowering | Yes (ir.Map with key, MappedStep with group_by) | implicit |
| `map_by` CWL transpilation | No (explicitly rejected) | `test_map_by_transpile_reports_explicit_grouping_not_supported` |
| `map_by` Nextflow/WDL transpilation | No (no Nextflow or WDL transpiler exists) | N/A |
| Table update semantics (`t // r`, `r // t`) | No (raises at forcing time) | not tested |

---

## 6. String Templates: Interpolation

Interpolation values represent string templates with variable substitution. They appear in output parameter defaults (`${outbase}.bam`), run parameter defaults, and potentially in step bindings.

### 6.1 Interpolation AST

```json
{
  "kind": "word",
  "parts": [
    { "kind": "literal", "text": "path/to/" },
    { "kind": "var", "name": "sample_id" },
    { "kind": "literal", "text": ".bam" }
  ]
}
```

Three part kinds:

| Kind | Fields | Semantics |
|------|--------|-----------|
| `literal` | `text` | A literal string segment, no substitution |
| `var` | `name` | A variable reference. The executor must resolve against the step's runtime binding context (workflow inputs + upstream step outputs). |
| `expr` | `text` | An arbitrary expression string (e.g., `${cpu * 2}`, `${outbase / 2}`). Evaluation is deferred. Executors may evaluate via shell, report an error, or reject the workflow. |

### 6.2 Resolution Rules

1. `literal` parts are concatenated verbatim.
2. `var` parts are resolved by looking up `name` in the runtime bindings dictionary (which maps parameter names to runtime values).
3. `expr` parts are target-dependent:
   - CWL: rejected at transpile time (unless the expression can be represented as a CWL parameter reference)
   - Shell executor: evaluated by the shell
   - Nextflow/WDL: evaluated via the target language's expression system if representable

### 6.3 Feature: Interpolation Implementation Status

| Sub-feature | Implemented | Tested |
|-------------|-------------|--------|
| Literal part | Yes | `test_output_glob_uses_cwl_expression` |
| Var part resolved in CWL globs | Yes (via `$(inputs.var + '.ext')`) | `test_output_glob_uses_cwl_expression` |
| Var part in run parameter defaults | Yes (serialized, resolved at compile time for built-in run params) | implicit |
| Expr part | Yes (serialized) but rejected by CWL transpiler | `test_rejects_output_expr_interpolation` |

---

## 7. Self-Containment Guarantee

The following information is embedded directly in the DAG JSON. No `path` field is ever used to re-read a source file.

| Information | Where in DAG | Example |
|-------------|-------------|---------|
| Bash script body | `steps[].script` | `"bwa mem -t ${cpu} ..."` |
| Input parameter specs | `steps[].inputs` | `{"fastq1": {"type": "file", "desc": "read 1"}}` |
| Output parameter specs | `steps[].outputs` | `{"bam": {"type": "file", "default": {"kind": "word", ...}}}` |
| Runtime resources | `steps[].run` | `{"cpu": {"type": "int", "value": 4}}` |
| Sub-workflow DAG | `steps[].definition.dag` (recursive) | full nested DAG |
| Task documentation | `steps[].task.doc` (embedded in task dict) | `"Align reads"` |
| Source file path (provenance only) | `steps[].path` | `"align.sh"` — never read |

---

## 8. Cross-Target Mapping Reference

| DAG concept | CWL | Nextflow | WDL |
|-------------|-----|----------|-----|
| Task step | `CommandLineTool` + `InitialWorkDirRequirement` | `process` with `script` block | `task` with `command {}` |
| Workflow step | `Workflow` + `steps[].run` → embedded `Workflow` | `workflow` + process `call` | `workflow` + `call` |
| Input binding | `source: "#main/name"` | `Channel.fromPath` / `take` | input variable |
| Step output binding | `source: "#main/step_id/out"` | process output channel | `call.step_name.output_name` |
| Literal binding | `default` on input port | `set` / `value` | `default` |
| Record merge | Flatten to individual field bindings | Flatten to tuple elements | Flatten to object fields |
| Field projection | `source: "#main/input_name/field"` | Channel `.map{ it.field }` | `input_name.field` |
| `map` (scatter) | `scatter` + `scatterMethod: dotproduct` | `each` + `fromFilePairs` | `scatter` block |
| `map_by` | Not supported | `groupTuple` + `splitText` | Manual grouping within scatter |
| CPU resource | `ResourceRequirement.coresMin` | `cpus` directive | `runtime { cpu: ... }` |
| Memory resource | `ResourceRequirement.ramMin` | `memory` directive | `runtime { memory: ... }` |
| Docker image | `DockerRequirement.dockerPull` | `container` directive | `runtime { docker: ... }` |
| Time resource | Not directly mappable | `time` directive | Not directly mappable |
| Interpolation (var) | `$(inputs.var + 'suffix')` | String interpolation (`"${params.var}"`) | String interpolation (`"${var}"`) |
| Interpolation (expr) | Rejected | Shell expression pass-through | Not supported |

---

## 9. Implementation Status Summary

### 9.1 Completed (no gaps)

- Task step emission with embedded script
- Workflow step emission with recursive embedded DAG
- Mapped task (`map f xs`) with input_schema / output_schema
- Mapped workflow (`map f xs` where f imports a `.swl`)
- Generated sub-workflow for `map` + lambda
- Root-level partial application of `map`
- Input, step_output, literal binding serialization
- Partial application with concrete record bindings
- CWL CommandLineTool emission (baseCommand, InitialWorkDirRequirement, inputs, outputs)
- CWL ResourceRequirement (cpu, memory)
- CWL DockerRequirement (image)
- CWL scatter with dotproduct scatterMethod
- CWL output glob from interpolation vars
- Interpolation AST serialization (literal, var)
- DAG round-trip (to_dict / from_dict)

### 9.2 Partially complete

| Feature | What works | What's missing |
|---------|------------|----------------|
| Input type/desc in DAG | Types merged from step specs if name matches | Inputs without matching step input names (lambda params, record fields) lack type/desc |
| Merge flattening in outputs | Chain desugaring produces flat outputs | Hand-written `//` trees not fully flattened during forcing |
| Merge binding serialization | Serialized as `{source: "merge"}` in DAG | No transpiler supports merge bindings (should be flattened before reaching transpilation) |
| Record binding serialization | Serialized as `{source: "record"}` in DAG | Records that don't immediately saturate a task reach the DAG; no transpiler supports them |
| Run parameter validation | Parameters with bound values are resolved | Parameters with no default and no bound value are silently dropped |
| Nested field projections | `Field(StepCall, ...)` and `Field(Input, ...)` in bindings | `Field(Field(...), ...)` nested chains not transpilable |
| Workflow definition embedding | Embedded as `definition.dag` | Top-level inputs/outputs/run fields duplicate info inside definition |

### 9.3 Not implemented

| Feature | Rationale |
|---------|-----------|
| DAG circularity validation | Missing validation step; should be added to `DAG.validate()` |
| `map_by` group_by transpilation | CWL rejected (no native grouped scatter); Nextflow/WDL not yet implemented |
| Table update semantics (`t // r`) | Rarely used; implementation deferred |
| `expr` interpolation in targets other than CWL | No Nextflow/WDL transpiler exists yet |
| Optional type (`?`) serialization | Optionality affects validation, not DAG structure; not yet needed by transpilers |
| Runtime/post-run-time checks (sec 6, 7) | Executor-level, not a DAG concern |
| Pre-run-time input validation wiring | Validator exists (`semantic/wf/validate.py`) but not connected to any execution path |
| Nextflow transpiler | Not built |
| WDL transpiler | Not built |
