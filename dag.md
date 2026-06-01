# SWL Compiled DAG JSON Specification

The DAG JSON is the **sole** output of the SWL compiler. It is a self-contained, minimal representation of a workflow as a directed acyclic graph of step invocations. Transpilers and execution engines must work entirely from DAG JSONs — no source files (`.swl`, `.sh`) may be read at transpilation or execution time.

---

## Design Principles

1. **Self-contained**: All task definitions, tool metadata, and scripts are embedded. A DAG JSON can be transpiled or executed without access to the original source files.
2. **Transpiler-neutral**: The schema represents function application, record merging, field projection, and batch semantics in a target-agnostic way.
3. **Minimal but sufficient**: Every field exists to serve at least one transpilation target or executor. Target-specific enrichment (e.g., CWL `scatterMethod`, Nextflow `fromChannel`) is added by the transpiler, never by the schema.

---

## Requirements due to Transpilation Targets

The DAG JSON is consumed by three transpiler targets (CWL, WDL 1.1, Nextflow DSL2), each with fundamentally different native primitives. This creates two hard requirements that the DAG schema and the compiler's IR-forcing phase must satisfy.

### R1. Explicit over inferred

Information that multiple transpilers independently re-derive must be made explicit in the DAG. Inference logic duplicated across N transpilers is a maintenance hazard — it means N copies of the same fragile heuristic, N places where bugs hide, and N times the work when the inference rule changes.

**Concrete examples of inference that should be eliminated:**

| Information | Currently inferred by | Why inference is fragile |
|---|---|---|
| Workflow output types | CWL transpiler looks up the source step's output type in `step.task.outputs[name].type`. WDL and Nextflow would need to do the same. | If the output source is an `Input` binding, a `Literal`, a `Record`, or a `Merge`, each transpiler must handle the type differently — or guess. Outputs without a direct step reference have no type at all. |
| Scatter vs broadcast ports in mapped steps | CWL transpiler checks whether a binding name appears in `map.source.columns`. WDL and Nextflow plans do the same check. | Requires understanding the internal structure of `map.source` and cross-referencing against `bindings`. If `map.source` is an `input` rather than a `table`, the port naming changes. Every transpiler must replicate this resolution logic. |
| Optional vs required parameters | Not serialized. CWL and WDL would each need to preserve the `?` from task annotation parsing through the entire compiler pipeline. | Currently lost at DAG serialization — the `?` flag is used during semantic checking but never written to the DAG JSON. Transpilers cannot recover it. |
| Tool identity for deduplication | CWL transpiler uses `step.id` as both call identity and tool identity. If two steps use the same underlying tool, they get separate `CommandLineTool` nodes. | CWL `$graph` requires one `CommandLineTool` per unique tool definition. Without an explicit `tool_id` or `definition_id`, the transpiler has no way to know that two steps share a definition. Currently this works because each step happens to have a unique task import, but it breaks when the same task is imported twice. |

**Rule:** If a piece of metadata would require code in every transpiler to recompute, it belongs in the DAG. The compiler has all the information — it should expose it once.

### R2. Flatten before emitting

Binding forms that no transpiler can handle natively must be resolved or flattened during IR forcing, not passed through to the DAG. The DAG is a contract: every binding form that appears in the output must be supportable by at least one target.

**Rationale:** Consider what happens when a `Merge` binding reaches the DAG:
- CWL transpiler: raises `ValueError` ("merge values are not supported")
- WDL transpiler: raises `ValueError` ("merge bindings must be flattened before transpilation")
- Nextflow transpiler: would attempt `.join()` on channels, but this only works when the merge represents a specific pattern (pipeline accumulation). For arbitrary merge trees the behavior is undefined.

The only sensible fix is to flatten merges in one place — the compiler — where all the semantic context is available. The transpiler is too late: by the time it sees the merge, the information about which fields overlap, which side should win, and whether the result is going to a single port has been lost.

**Binding forms covered by this requirement:**

| Binding form | Transpiler support | Required action |
|---|---|---|
| `merge` (`r1 // r2`) | Rejected by CWL and WDL. Nextflow `.join()` is fragile and pattern-specific. | Flatten into per-field bindings during forcing. If flattening is impossible, raise compile-time error. |
| `record` (non-saturating) | CWL rejects. WDL and Nextflow can emit struct/tuple construction, but the transpiler must know the record's shape. | **Preferred:** saturate into task input bindings during forcing (each field becomes a separate step input). **Fallback:** emit with enough metadata for transpiler-side construction. |
| `Field(Field(...))` (nested) | CWL rejects. WDL and Nextflow can technically resolve the chain, but the binding tree is complex. | **Preferred:** resolve to a direct step-output reference if the compiler knows the intermediate type. **Fallback:** Emit with a chain length hint so the transpiler can generate resolution code. |
| `table` with complex column bindings | CWL supports it. WDL and Nextflow need to reconstruct arrays from step outputs, not just workflow inputs. | The table source should resolve to concrete array bindings during forcing. Column bindings that reference `Field(MappedStep, ...)` should be traceable to a single array-valued source. |

**Rule:** The forcing phase is the last point where the compiler has complete type information and binding context. It should use that context to resolve, flatten, or reject binding forms that would require context-free transpilers to guess. Every merge or record that reaches the DAG is a bug in the compiler, not a limitation of the transpiler.

---

## Implementation priorities

| Priority | What to do | Affected files |
|---|---|---|
| P0 | Add merge flattening pass to `force.py` | `ir/force.py`, `ir/dag.py` |
| P1 | Serialize `map.scatter` / `map.broadcast` in `MappedStep` | `ir/dag.py`, `ir/force.py` |
| P1 | Serialize `optional` flag in `InputSpec` / `ParamSpec` | `ir/dag.py` |
| P2 | Serialize output `type` in top-level `outputs` | `ir/dag.py`, update transpilers to read typed format |
| P2 | Flatten non-saturating records or add shape metadata | `ir/force.py` |
| P3 | Add `tool_id` / `definition_id` for deduplication | `ir/dag.py`, transpilers |

---

## 1. Top-Level Structure

```json
{
  "inputs": { "<name>": <InputSpec>, ... },
  "steps": [ <StepSpec>, ... ],
  "outputs": { "<name>": { "type": <type string>, "desc": <string or null>, "value": <Binding> }, ... }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `inputs` | object | yes | Map of workflow input name → InputSpec |
| `steps` | array | yes | Ordered list of step specifications |
| `outputs` | object | yes | Map of workflow output name to output descriptor |

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

Optionality is serialized via the `"optional"` field on `InputSpec` and `ParamSpec` (§2.2, §2.3). When absent, the field is treated as required.

### 2.2 InputSpec

```json
{
  "type": "file" | "str" | "int" | "float" | "[file]" | "[str]" | "[int]" | "[float]" | null,
  "desc": "description string or null",
  "optional": true
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `type` | string or null | yes | Scalar or array type |
| `desc` | string or null | no | Human-readable description |
| `optional` | boolean | no | If true, the input may be omitted at runtime. Default: false. |

Declares a named workflow input.

### 2.3 ParamSpec

Describes an input or output parameter of a task or workflow step.

```json
{
  "type": "file" | "str" | "int" | "float" | "[file]" | "[str]" | "[int]" | "[float]" | null,
  "desc": "human-readable description or null",
  "default": <Interpolation> or null,
  "optional": true
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `type` | string or null | yes | Scalar or array type |
| `desc` | string or null | no | Human-readable description |
| `default` | Interpolation or null | no | Default value expression. For output parameters, this is the expected file path or glob pattern. |
| `optional` | boolean | no | If true, the parameter may be omitted. Default: false. |

**Transpiler requirements:**
- CWL: optional types emit as `["null", <type>]` (union with null)
- WDL: optional types emit with `?` suffix (e.g., `File?`)
- Nextflow: optional inputs use `optional` qualifier or `?:` default expression

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
Used by pipeline desugaring (`A | B | C`) to accumulate outputs.

**Contract: Merge bindings must not appear in the final DAG.** The compiler's IR forcing phase must flatten merge trees into flat per-field bindings before emitting the DAG. This is because no transpiler target (CWL, WDL, Nextflow) has a native record-merge primitive — each would need to implement decomposition differently. The flattening is done once in the compiler, not N times in N transpilers.

If a merge binding cannot be flattened (e.g., merging two non-overlapping record-valued step outputs into a single port), the compiler must raise a compile-time error rather than emitting an unsupported binding form. See §4.2 for current status.

**`record`**: A constructed record literal.
```json
{ "source": "record", "fields": { "<name>": <Binding>, ... } }
```
Used when partial application creates a concrete input record or when a lambda produces a result record.

**Contract: Record bindings should saturate a task call during forcing.** A record that directly feeds a task's input ports should be flattened into individual bindings (each field becomes a separate step input). Non-saturating records (those that flow through intermediate bindings without being consumed by a concrete tool) reach the DAG and require transpiler-side struct construction. Transpiler guidance:
- CWL: emit an `ExpressionTool` that constructs the record from individual inputs
- WDL: emit a `struct` definition + struct literal
- Nextflow: emit a `tuple` or channel-of-records

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
| `map.scatter` | array of strings | The subset of input parameter names that are scatter ports — their values come from table columns and vary per row. Must be a subset of `bindings` keys. |
| `map.broadcast` | array of strings | The subset of input parameter names that are broadcast ports — their values are constant across all rows. Must be a subset of `bindings` keys. |

The `scatter` and `broadcast` fields disambiguate how each binding is applied during batch execution. Every input parameter must appear in exactly one of `scatter` or `broadcast`. If absent, transpilers must infer the split by analyzing `map.source` column bindings — this inference is fragile and duplicates compiler logic. The compiler should always populate both fields.

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
    "dag": { "inputs": {...}, "steps": [...], "outputs": {...} },
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
| Mapped step with `map_by` (group_by) | Yes (compiler produces group_by) | implicit |
| `map.scatter` / `map.broadcast` port disambiguation | No (must be inferred by transpilers) | — |
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

Same structure as `map` but with `map.group_by` set to the grouping key column name.

**Transpiler contract:** Each transpiler target has a different native mechanism for grouped scatter:

| Target | Mechanism | Plan reference |
|--------|-----------|----------------|
| CWL | Synthetic `ExpressionTool` (grouping) + scattered wrapper over group JSONs | `map_by_cwl.md` |
| WDL 1.1 | `collect_by_key()` function from standard library | `wdl.md §5` |
| Nextflow | Channel `.groupTuple()` operator | `nf.md §5` |

The DAG's `map.group_by` field is the sole signal. Each transpiler maps it to the target's native grouping primitive. No transpiler should implement grouping logic itself — they generate code that delegates to the target runtime.

### 5.3 Feature: Batch Semantics Implementation Status

| Sub-feature | Implemented | Tested |
|-------------|-------------|--------|
| `map` on task | Yes | `test_batch_mapped_task_emits_scatter_and_tab_column_input_type` |
| `map` on workflow | Yes | `test_batch_mapped_workflow_emits_scattered_subworkflow` |
| `map` on lambda (generated sub-workflow) | Yes | `test_batch_mapped_simple_lambda_emits_generated_scattered_subworkflow` |
| `map` on partial application (root `map f`) | Yes | `test_root_partial_map_transpiles_as_scattered_subworkflow` |
| `map_by` compile-time lowering | Yes (ir.Map with key, MappedStep with group_by) | implicit |
| `map.scatter` / `map.broadcast` port lists | No (transpilers infer split; `map_by` plan may require this) | — |
| `map_by` CWL transpilation | No (explicitly rejected — planned via ExpressionTool + scatter) | `test_map_by_transpile_reports_explicit_grouping_not_supported` |
| `map_by` WDL transpilation | No (planned via `collect_by_key()`) | Not built |
| `map_by` Nextflow transpilation | No (planned via `groupTuple()`) | Not built |
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

| Target | `expr` handling |
|--------|-----------------|
| CWL | Rejected by default. If `InlineJavascriptRequirement` is acceptable, emit as `$(<text>)` passthrough. |
| Shell executor | Evaluated by the shell (the expression text is valid bash). |
| WDL 1.1 | Emit as `~{<text>}` — WDL 1.1 expression syntax supports most arithmetic and string operations. |
| Nextflow | Emit as `${ <text> }` — Groovy expression passthrough. |

**Transpiler flag for `expr`:** The DAG may contain `expr` parts that reference variables outside the interpolation context (e.g., `${cpu * 2}` references `cpu`, which is a run parameter, not a binding). Transpilers should validate variable availability within their target's scope before emitting `expr` passthrough. If a variable referenced in `expr` is not available in the target's expression scope, the transpiler must reject with a clear error.

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

| DAG concept | CWL | Nextflow | WDL 1.1 |
|-------------|-----|----------|---------|
| Task step | `CommandLineTool` + `InitialWorkDirRequirement` | `process` with `script` block | `task` with `command <<< >>>` |
| Workflow step | `Workflow` + `steps[].run` → embedded `Workflow` | `workflow` + process `call` | `workflow` + `call` |
| Input binding | `source: "#main/name"` | `Channel.fromPath` / `take` | input variable |
| Step output binding | `source: "#main/step_id/out"` | process output channel | `call.step_name.output_name` |
| Literal binding | `default` on input port | `Channel.value(val)` | inline literal |
| Record merge | Flatten at compiler time — not supported in DAG | Flatten at compiler time — not supported in DAG | Flatten at compiler time — not supported in DAG |
| Record binding | `ExpressionTool` construction | `tuple` / `record` construction | `struct` literal |
| Field projection | `source: "#main/input_name/field"` | Channel `.map{ it.field }` | `call.output.field` |
| `map` (scatter) | `scatter` + `scatterMethod: dotproduct` | Implicit scatter via channel `.join()` | `scatter` block |
| `map_by` (planned) | `ExpressionTool` grouping + scattered wrapper | `.groupTuple()` | `collect_by_key()` + scatter |
| CPU resource | `ResourceRequirement.coresMin` | `cpus` directive | `cpu:` in `requirements {}` |
| Memory resource | `ResourceRequirement.ramMin` | `memory` directive | `memory:` in `requirements {}` |
| Docker image | `DockerRequirement.dockerPull` | `container` directive | `container:` in `requirements {}` |
| Time resource | Hint (non-standard) | `time` directive | `time_minutes:` in `requirements {}` (WDL 1.1) |
| Optional type (`file?`) | `["null", "File"]` | `optional` qualifier | `File?` suffix |
| Interpolation (var) | `$(inputs.var + 'suffix')` | `"${params.var}"` | `~{var}` |
| Interpolation (expr) | Rejected or `InlineJavascriptRequirement` | Shell expression pass-through | `~{expr}` passthrough |

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
| Merge flattening in outputs and bindings | Chain desugaring produces flat outputs | Hand-written `//` trees not fully flattened during forcing; merge bindings still reach DAG |
| Record binding handling | Serialized as `{source: "record"}` in DAG | Non-saturating records reach DAG; CWL rejects, WDL/Nextflow need struct/record construction |
| Run parameter validation | Parameters with bound values are resolved | Parameters with no default and no bound value are silently dropped |
| `map.scatter` / `map.broadcast` port lists | Not serialized | Transpilers must infer scatter vs broadcast split from `map.source` bindings — fragile |
| Nested field projections | `Field(StepCall, ...)` and `Field(Input, ...)` in bindings | `Field(Field(...), ...)` nested chains not transpilable |
| Workflow definition embedding | Embedded as `definition.dag` | Top-level inputs/outputs/run fields duplicate info inside definition |
| Output type serialization | Version `"1.0"`: no output types | Version `"1.1"`: outputs include explicit `type` (proposed, not implemented) |
| Optional type serialization | Not serialized | Proposed: `"optional": true` on InputSpec/ParamSpec |

### 9.3 Not implemented

| Feature | Rationale |
|---------|-----------|
| DAG circularity validation | Missing validation step; should be added to `DAG.validate()` |
| Merge flattening guarantee | Force phase does not yet flatten all merge bindings; transpilers still encounter them |
| `map_by` transpilation (CWL) | Planned: ExpressionTool grouping + scattered wrapper (`map_by_cwl.md`) |
| `map_by` transpilation (WDL) | Implemented: `collect_by_key()` + scatter (`wdl.md §5`) |
| `map_by` transpilation (Nextflow) | Implemented: `groupTuple()` (`nf.md §5`) |
| `map.scatter` / `map.broadcast` serialization | Required by all three transpiler plans; not yet emitted by compiler |
| Record binding transpilation | CWL: ExpressionTool; WDL: struct literal; Nextflow: tuple/channel construction |
| Expr interpolation (`expr` kind) in CWL | Planned: optional `InlineJavascriptRequirement` passthrough |
| Time resource in CWL | Planned: emit as hints entry |
| Table update semantics (`t // r`) | Rarely used; implementation deferred |
| Nextflow transpiler | Implemented (`nf.md`, `python/swl/transpile/nf/`) — Phases 1–8 complete |
| WDL transpiler | Implemented (`wdl.md`, `python/swl/transpile/wdl/`) — Phases 1–7 complete |
