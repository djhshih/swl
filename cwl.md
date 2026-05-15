# Plan: transpile compiled SWL DAG JSON to packed CWL

This note describes:
- how to map the current compiled SWL JSON into packed CWL
- what additional information CWL generation needs
- what should change in the compiled SWL JSON so CWL transpilation is straightforward and reliable

Target:
- emit a single packed CWL document with `cwlVersion` and `$graph`
- include one `Workflow` node for the compiled SWL workflow
- include one `CommandLineTool` node per distinct task definition used by the workflow
- use `InitialWorkDirRequirement` to materialize each task script into the execution directory

Reference style examined:
- `cipher.pack.cwl`

---

## 1. What the current compiled SWL JSON already gives us

The current DAG JSON already contains several things needed for CWL generation:

- top-level workflow inputs
- top-level workflow outputs
- task call instances with dependency order
- task input bindings
- task output declarations
- normalized run values
- embedded task script text
- task source path and task name

This is enough to generate a first packed CWL file for simple workflows.

In particular, each compiled task already contains enough information to generate a `CommandLineTool` body:
- interface
- script body
- resource-like run fields (`cpu`, `memory`, `time`, `image`)

And each task call already contains enough information to generate a workflow step:
- step id
- step inputs
- step outputs
- edges to upstream steps or workflow inputs

---

## 2. Packed CWL structure to generate

The generated CWL should look like:

```json
{
  "cwlVersion": "v1.0",
  "$graph": [
    { "id": "#align", "class": "CommandLineTool", ... },
    { "id": "#sort", "class": "CommandLineTool", ... },
    { "id": "#function", "class": "Workflow", ... }
  ]
}
```

Target version:
- `v1.0`

The transpiler should emit CWL `v1.0` to match the style of `cipher.pack.cwl`.

The packed file should contain:

1. one `CommandLineTool` entry per distinct task definition
2. one `Workflow` entry representing the compiled DAG

Distinct task definitions should be deduplicated by task path or another stable task identity, not by step id.

---

## 3. Mapping SWL task definitions to CWL CommandLineTool

## 3.1 Tool identity

For each imported task definition, create a packed tool entry.

Tool id rule:
- derive the tool id from the variable name bound to the task import in the SWL source
- example:

```swl
align = import "align.sh"
sort  = import "sort.sh"
```

should produce packed tool ids:
- `#align`
- `#sort`

This matches the language-level naming the user wrote, and gives stable readable ids.

Implication:
- the compiled artifact must preserve the import binding name for each imported task definition
- the transpiler should not invent tool ids from path basenames unless no import-bound name is available

---

## 3.2 Base command

Since SWL tasks currently embed shell script bodies, the simplest CWL lowering is:

```json
"baseCommand": ["bash", "script.sh"]
```

and then use:

```json
{
  "class": "InitialWorkDirRequirement",
  "listing": [
    {
      "entryname": "script.sh",
      "entry": "...task shell body..."
    }
  ]
}
```

This matches the style in `cipher.pack.cwl`.

Notes:
- later we may want `ShellCommandRequirement` and inline command fragments, but script materialization is simpler and closer to current SWL task semantics
- one script filename per tool is enough if the file is local to the tool execution directory

---

## 3.3 Tool inputs

Each compiled task input should become a CWL input parameter.

Example current SWL task metadata:

```json
"inputs": {
  "bam": {"type": "file", "desc": "input bam"},
  "outbase": {"type": "str", "desc": "output base name"}
}
```

Possible CWL lowering:

```json
"inputs": [
  {"id": "#tool/sort/bam", "type": "File", "doc": "input bam"},
  {"id": "#tool/sort/outbase", "type": "string", "doc": "output base name"}
]
```

Type mapping plan:

- `file` -> `File`
- `str` -> `string`
- `int` -> `int`
- `float` -> `float` if later added
- `bool` -> `boolean` if later added
- `memory` -> likely `int` in MiB for resource fields, not a normal CLI file/string input
- `time` -> likely `int` in minutes for resource fields, not a normal CLI file/string input

For now, task `run` fields should not necessarily become ordinary CWL formal inputs unless they are intended to be overridable at workflow runtime.

---

## 3.4 Tool outputs

Each compiled task output should become a CWL output parameter.

Current SWL compiled JSON already includes output defaults in interpolation form, for example:

```json
"outputs": {
  "bam": {
    "type": "file",
    "default": {
      "kind": "word",
      "parts": [
        {"kind": "var", "name": "outbase"},
        {"kind": "literal", "text": ".bam"}
      ]
    }
  }
}
```

This can lower to:

```json
{
  "id": "#tool/sort/bam",
  "type": "File",
  "outputBinding": {
    "glob": "$(inputs.outbase + '.bam')"
  }
}
```

So the transpiler needs a conversion from SWL interpolation JSON to CWL expression strings.

That is feasible for:
- literal pieces
- variable references
- simple concatenation

Potentially not feasible yet for arbitrary `${expr}` fragments unless we define a translation subset.

---

## 3.5 Tool resource requirements

Current SWL run fields:
- `cpu`
- `memory`
- `time`
- `image`

Reasonable CWL lowering:

### cpu

```json
{
  "class": "ResourceRequirement",
  "coresMin": 4
}
```

### memory

Current SWL semantic normalization is MiB.
That maps well to:

```json
{
  "class": "ResourceRequirement",
  "ramMin": 8192
}
```

### time

CWL core does not have a standard runtime time limit field.
Options:
- ignore in first pass
- preserve as extension metadata
- emit in `hints`

Recommended first-pass behavior:
- preserve `time` in a non-standard `hints` or extension field
- do not pretend it is standard CWL runtime semantics

### image

This wants container lowering:

```json
{
  "class": "DockerRequirement",
  "dockerPull": "djhshih/seqkit:0.1"
}
```

So SWL `run.image` can map cleanly to CWL `DockerRequirement`.

---

## 4. Mapping SWL DAG task calls to CWL Workflow steps

Each compiled `tasks[]` item becomes one CWL workflow step.

Example SWL step-like data:

```json
{
  "id": "t2",
  "name": "sort",
  "inputs": {
    "bam": {"source": "task", "task": "t1", "output": "bam"},
    "outbase": {"source": "input", "name": "outbase"}
  },
  "outputs": {...}
}
```

Lower to CWL step:

```json
{
  "id": "#main/t2",
  "run": "#tool/sort",
  "in": [
    {"id": "#main/t2/bam", "source": "#main/t1/bam"},
    {"id": "#main/t2/outbase", "source": "#main/outbase"}
  ],
  "out": ["#main/t2/bam", "#main/t2/bai"]
}
```

The `deps` field is useful for validation and ordering, but not necessary for final CWL edges because CWL step dependencies are derived from `source` links.

---

## 5. Mapping top-level SWL DAG inputs to CWL Workflow inputs

Current compiled JSON top-level inputs:

```json
"inputs": {
  "fastq1": {"type": "file", "desc": "paired-end reads"},
  "outbase": {"type": "str", "desc": "output base name"}
}
```

Lower to workflow inputs:

```json
"inputs": [
  {"id": "#main/fastq1", "type": "File", "doc": "paired-end reads"},
  {"id": "#main/outbase", "type": "string", "doc": "output base name"}
]
```

This part is already fairly direct.

---

## 6. Mapping top-level SWL DAG outputs to CWL Workflow outputs

Current compiled JSON top-level outputs reference:
- task outputs
- sometimes input/literal/record/merge/function-shaped values

For CWL workflow outputs, the easiest supported case is:
- output bound directly to a step output

Example:

```json
"outputs": {
  "bam": {"source": "task", "task": "t2", "output": "bam"}
}
```

Lower to:

```json
{
  "id": "#main/bam",
  "type": "File",
  "outputSource": "#main/t2/bam"
}
```

Harder cases:
- record-valued outputs
- merged outputs
- function-valued outputs
- literal outputs that are not files

For first-pass CWL transpilation, we should restrict support to compiled DAGs whose final outputs are flat named bindings to:
- workflow inputs
- task outputs
- maybe literals

If compiled output still contains `source: function`, CWL transpilation should fail clearly.

---

## 7. InitialWorkDirRequirement usage

Use `InitialWorkDirRequirement` in each generated `CommandLineTool`.

Planned form:

```json
{
  "class": "InitialWorkDirRequirement",
  "listing": [
    {
      "entryname": "script.sh",
      "entry": "...shell script text..."
    }
  ]
}
```

This is appropriate because SWL compiled JSON already embeds the exact shell body needed for execution.

Longer-term option:
- also materialize small helper files if SWL later grows that concept

---

## 8. What in the current SWL JSON is awkward for CWL transpilation

The current compiled JSON is close, but some parts are still too executor-internal or too syntax-shaped.

Main issues:

1. task outputs store SWL interpolation syntax trees, not already-lowered filename expressions
2. task inputs do not distinguish ordinary CLI inputs from run/resource controls
3. top-level outputs can still contain non-CWL-friendly binding forms
4. task command construction relies on implicit shell variable naming rather than explicit input binding metadata
5. there is no explicit per-task CWL-friendly command template
6. file-vs-value distinctions are present, but not rich enough for all CWL lowering decisions

---

## 9. Recommended changes to compiled SWL JSON

This section is the main design recommendation.

## 9.1 Add explicit task interface sections by role

Current task entry shape roughly has:
- `inputs`
- `outputs`
- `run`
- `script`

Keep that, but make each parameter more explicit about execution role.

Recommended normalized param schema:

```json
{
  "type": "file",
  "desc": "input bam",
  "role": "input"
}
```

and for outputs:

```json
{
  "type": "file",
  "desc": "output alignment",
  "role": "output",
  "path": {
    "engine": "swl_interp",
    "value": {...}
  }
}
```

and for run params:

```json
{
  "type": "int",
  "desc": null,
  "role": "resource",
  "value": 4,
  "resource": "cpu"
}
```

Why:
- CWL lowering needs to know whether a field is a real workflow/CLI input, a declared output path, or a resource/container setting
- right now this is implied by section location, but making it explicit simplifies a separate transpiler

---

## 9.2 Add explicit output path expressions

Current output `default` values are SWL interpolation trees.
That is good raw information, but the transpiler should not have to infer semantics from generic `default`.

Recommended change:
- rename task output `default` to something output-specific, such as `path` or `glob`

Example:

```json
"outputs": {
  "bam": {
    "type": "file",
    "desc": "output alignment",
    "path": {
      "kind": "word",
      "parts": [
        {"kind": "var", "name": "outbase"},
        {"kind": "literal", "text": ".bam"}
      ]
    }
  }
}
```

Better yet, include both preserved SWL form and pre-rendered target-friendly forms:

```json
"outputs": {
  "bam": {
    "type": "file",
    "desc": "output alignment",
    "path": {
      "swl": {...},
      "cwl": "$(inputs.outbase + '.bam')"
    }
  }
}
```

Why:
- CWL `outputBinding.glob` is path-oriented
- task output defaults in SWL are not really defaults in the same sense as input defaults
- explicit path semantics make the contract clearer

---

## 9.3 Add explicit task command description, not just raw script

Current compiled task JSON stores only:

```json
"script": "...bash text..."
```

For CWL, that is enough if we always use `bash script.sh`.
But it would help to carry an explicit execution stanza:

```json
"execution": {
  "kind": "bash_script",
  "script_name": "script.sh",
  "script": "...bash text..."
}
```

Why:
- makes the compiled JSON less ad hoc
- separates execution model from arbitrary top-level field naming
- gives future room for non-bash task kinds

---

## 9.4 Add explicit CWL-friendly value kinds for bindings

Current binding forms are:
- `input`
- `task`
- `literal`
- `record`
- `merge`
- `field`
- `task_call`
- `function`

For CWL transpilation, the supported subset should be made explicit.

Recommended addition:
- define and document which output/input binding forms are valid for CWL lowering
- optionally add a derived normalized subset, e.g.:

```json
"cwl_binding": {
  "kind": "step_output",
  "step": "t2",
  "output": "bam"
}
```

or

```json
"cwl_binding": {
  "kind": "workflow_input",
  "name": "outbase"
}
```

Why:
- a transpiler should not need to reverse-engineer arbitrary merge trees if the DAG finalization phase can already flatten them
- forcing is the best place to reject or normalize non-CWL-friendly output shapes

---

## 9.5 Add explicit output types at top-level workflow outputs

Current top-level outputs are only bindings.
For CWL workflow outputs, a type is required.

Recommended compiled output shape:

```json
"outputs": {
  "bam": {
    "type": "file",
    "desc": "output alignment",
    "value": {"source": "task", "task": "t2", "output": "bam"}
  }
}
```

instead of only:

```json
"outputs": {
  "bam": {"source": "task", "task": "t2", "output": "bam"}
}
```

Why:
- CWL workflow outputs require declared types
- inferring them again during transpilation is possible, but unnecessary duplication

This is one of the most important JSON changes.

---

## 9.6 Distinguish task-definition identity from task-call identity

Current JSON has:
- task call `id` like `t1`
- task `name`
- task `path`

For packed CWL we need two identities:
- definition id for the deduplicated `CommandLineTool`
- call id for each workflow step instance

Recommended task call addition:

```json
{
  "id": "t2",
  "tool": "tool/sort",
  "name": "sort",
  ...
}
```

and top-level packed-tool table, or enough stable identity to derive it cleanly.

Why:
- several steps may refer to the same tool definition
- packed CWL needs a separate reusable tool id

---

## 9.7 Preserve container/resource info in a normalized role-aware form

Current `run` values are useful, but for CWL lowering we should make standard resource mappings explicit.

Recommended normalized shape:

```json
"run": {
  "cpu": {"type": "int", "value": 4, "cwl_resource": "coresMin"},
  "memory": {"type": "memory", "value": 8192, "unit": "MiB", "cwl_resource": "ramMin"},
  "image": {"type": "str", "value": "ubuntu:22.04", "cwl_requirement": "DockerRequirement"},
  "time": {"type": "time", "value": 30, "unit": "minutes", "cwl_hint": "timeLimit"}
}
```

Why:
- avoids hard-coding SWL built-in run names in the transpiler
- makes the JSON more target-aware without fully baking in CWL

---

## 9.8 Add a compiled-workflow name/id

Workflow id rule:
- derive the CWL workflow id from the file stem of the source `.swl` file
- example:
  - `tests/function.swl` -> `#function`
  - `tests/pipe.swl` -> `#pipe`

This should be the primary rule for the packed workflow id.

Recommended compiled top-level additions:

```json
{
  "source": {
    "path": "tests/function.swl",
    "stem": "function"
  },
  "name": "function",
  ...
}
```

Why:
- packed CWL needs a stable workflow id
- the chosen rule is now explicit: use the `.swl` file stem
- preserving the source stem in the compiled artifact keeps the transpiler deterministic and avoids recomputing naming policy elsewhere

---

## 10. Proposed transpilation algorithm

## Phase A: load and validate compiled DAG JSON

1. load JSON
2. validate supported binding/output forms for CWL lowering
3. ensure no top-level output is function-valued
4. ensure task definitions are present and self-contained

Reject clearly if unsupported:
- `source: function`
- complex record/merge outputs not flattened to named workflow outputs
- output path expressions outside supported interpolation subset

---

## Phase B: build packed tool table

1. collect distinct task definitions by stable task identity
2. for each unique task definition:
   - create `CommandLineTool`
   - map inputs
   - map outputs
   - create `InitialWorkDirRequirement`
   - map `run` values into `ResourceRequirement`, `DockerRequirement`, and optional hints

---

## Phase C: build workflow node

1. create workflow inputs from top-level DAG inputs
2. create a step for each compiled task call
3. wire each step input from:
   - workflow input
   - upstream step output
   - literal/default strategy if supported
4. create workflow outputs from top-level DAG outputs

---

## Phase D: emit packed CWL

Emit:

```json
{
  "cwlVersion": "v1.0",
  "$graph": [tools..., workflow]
}
```

Tool nodes may appear before or after workflow; either is fine as long as ids resolve.

---

## 11. Important first implementation constraints

To keep the first transpiler tractable, the first supported subset should be:

1. task bodies are shell scripts lowered through `InitialWorkDirRequirement`
2. workflow outputs must be flat named bindings
3. step input sources must be:
   - workflow input
   - upstream task output
   - maybe scalar literal
4. output path interpolation must stay within:
   - literals
   - `${name}` variable substitution
   - concatenation
5. only one top-level workflow is emitted
6. no attempt yet to lower lazy function-valued outputs to CWL

That subset already covers the current `pipe/function/explicit` style examples.

---

## 12. Recommended implementation steps

1. implement the transpiler under:
   - `python/swl/transpile/cwl/`

   Recommended layout:
   - `python/swl/transpile/cwl/__init__.py`
   - `python/swl/transpile/cwl/emit.py`
   - `python/swl/transpile/cwl/types.py`
   - `python/swl/transpile/cwl/cli.py` if a dedicated CLI is added

2. add a CLI:
   - `python -m swl.cwl tests/function.swl`
   - or transpile from compiled JSON directly

3. decide the entrypoint:
   - preferred: transpile from compiled DAG JSON, not directly from `.swl`
   - this keeps CWL generation downstream of forcing and avoids duplicating semantic/lowering logic

4. first make forcing emit CWL-friendly metadata additions:
   - workflow id/name
   - top-level output types
   - explicit output path field
   - stable tool identity

5. implement a strict transpiler for the supported subset
6. add golden tests comparing generated packed CWL JSON structure

---

## 13. Summary of required SWL JSON changes

Minimum strongly recommended changes:

1. **top-level workflow metadata**
   - add workflow source metadata including `.swl` path and file stem
   - derive the CWL workflow id from the `.swl` file stem

2. **top-level outputs**
   - change from raw binding only to:
     - `type`
     - optional `desc`
     - `value`

3. **task output metadata**
   - rename/clarify output `default` as output path semantics
   - ideally expose `path` instead of generic `default`

4. **task/tool identity**
   - preserve the variable name bound to each imported task
   - derive the packed `CommandLineTool` id from that import-bound variable name
   - keep per-call task ids separate from per-definition tool ids

5. **task execution metadata**
   - wrap `script` in an explicit execution object, or at least standardize script file metadata

6. **run/resource metadata**
   - preserve normalized values
   - optionally add explicit target-neutral mapping hints like `resource: cpu`, `unit: MiB`

7. **CWL-friendly output normalization**
   - ensure forcing can flatten final outputs into simple named values suitable for CWL outputSource lowering

These changes would make a packed-CWL transpiler much simpler and less inference-heavy.

---

## 14. Recommended design principle

Do not make the CWL transpiler recover hidden semantics from low-level executor JSON if forcing can expose them directly.

Best boundary:
- SWL semantic/lowering/forcing should produce a self-contained compiled artifact with explicit task interface, output-path, resource, and workflow-output metadata
- CWL transpilation should mostly be a deterministic format conversion on that compiled artifact

That keeps CWL support maintainable and avoids duplicating compiler logic in the transpiler.
