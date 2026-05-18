## Core decision

Tasks and workflows are record-oriented, but with batch support there are now distinct categories that should be described separately.

- `rec`: `{ str: file|num|str, ... }`
- `[rec]`: homogeneous array of records
- **task**: `rec -> rec`
- **simple workflow**: `rec -> rec`
- **batch workflow**: either `[rec] -> rec` or `[rec] -> [rec]`

For the first implementation, `map` creates batch behavior. A mapped callee may be:
- a task (`rec -> rec`)
- a simple workflow (`rec -> rec`)

For now, mapping over a batch workflow is forbidden.
That is, if a workflow already takes `[rec]`, then `map` may not be applied to it.

## 1. Task

A task is the existing shell-script unit described by task annotations.

- signature: `rec -> rec`
- compiled JSON representation: one task entry in `tasks[]`
- this does **not** change semantically for the new feature

Example task:

```sh
#@  Sort alignment by coordinates and index
# in
#   bam      file
#   outbase  str
# out
#   bam      file = ${outbase}.bam
#   bai      file = ${outbase}.bai

samtools sort ...
samtools index ...
```

## 2. Simple workflow

This is the existing SWL model.

Example:

`call_variant.wf`
```swl
align = import "align.sh"
sort = import "sort.sh"
call = import "call.sh"
align | sort | call
```

- signature: `rec -> rec`
- compiled JSON representation: exactly the current DAG form, i.e. a flat `tasks[]` list plus top-level `inputs` and `outputs`
- examples: `tests/dag/pipe.json`, `tests/dag/function.json`, `tests/dag/partial.json`
- this does **not** change semantically for the new feature

Current shape, as seen in `tests/dag/pipe.json`:
- each concrete task call appears in top-level `tasks[]`
- each task entry has `id`, `path`, `deps`, `inputs`, `bindings`, `outputs`, `run`, `script`
- top-level `outputs` point to concrete task outputs, e.g. `{ "source": "task", "task": "call", "output": "bcf" }`

So today, a simple workflow is compiled by fully materializing its internal task DAG into the top-level `tasks[]` array.

## 3. Batch workflow

A batch workflow is a workflow whose parameter is `[rec]`.
For the intended final target, a batch workflow may return either:
- `rec`
- `[rec]`

Example reduction case (`[rec] -> rec`):

```sh
#@  Merge mutations
# in
#   bcf      [file]
#   outbase  str
# out
#   bcf      file  =  ${outbase}.bcf

bcftools merge -m both -o ${outbase}.bcf ${bcf}
```

```panel.swl
call_variant = import "call_variant.wf"
merge = import "merge.sh"

\xs -> 
   calls = map call_variant xs
   merge { bcf: calls.bcf }
```

Example passthrough case (`[rec] -> [rec]`):

```panel_only.swl
call_variant = import "call_variant.wf"

\xs ->
   map call_variant xs
```

Remarks
- `map` lifts a function from `rec -> rec` to `[rec] -> [rec]`
- `panel.swl` is a batch workflow with signature `[rec] -> rec`
- `panel_only.swl` is a batch workflow with signature `[rec] -> [rec]`
- `calls` is a value of type `[rec]`
- get operator `.` extends to `[rec]` and returns an array: `calls.bcf` returns `[file]`
- because `map` applies one statically typed `rec -> rec` function uniformly to all elements, `[rec]` produced by `map` is homogeneous, not heterogeneous
- if `field` is not present in the mapped element record type, `xs.field` is a compile-time error
- unlike a simple workflow, a batch workflow should not always be compiled by expanding the mapped workflow into many concrete calls; instead, mapped execution must remain symbolic until target lowering so that runtime cardinality can be preserved

## New builtins and operators

### `map`

`map` is a builtin processed in the same builtin-handling path as `import`.
Syntactically it still appears as ordinary application, but semantically it is builtin.

`map` takes a function `f1` with signature `rec -> rec` and returns a function
`f2` with signature `[rec] -> [rec]`.

Allowed mapped callees in the target design:
- imported task
- imported simple workflow
- lambda of type `rec -> rec`
- partially applied function whose remaining signature is `rec -> rec`

Forbidden for now:
- `map` applied to a batch workflow (`[rec] -> rec` or `[rec] -> [rec]`)
- nested batch mapping where the mapped callee itself consumes `[rec]`

When `f2` is applied to a record array `[rec]`, `f1` is run on each record of the array,
producing output records that are then structured into a homogeneous record array.

### `.` on `[rec]`

If `xs : [rec]`, then:

```swl
xs.field
```

returns `[t]`, where `field : t` in each element record.

If `field` is absent from the statically known element record type, this is a compile-time error.

## Concrete implementation target

This section describes the final target behavior, not just an intermediate sketch.
The implementation may be phased, but the end state should support compilation and CWL transpilation of batch workflows into runnable CWL.

This target is based on the current compiler shape in `python/swl`:
- lowering builds IR around scalar/record workflows (`python/swl/ir/lower.py`, `python/swl/ir/node.py`)
- forcing materializes execution DAG JSON (`python/swl/ir/force.py`, `python/swl/ir/dag.py`)
- CWL transpilation currently assumes only simple workflows and flat task-output wiring (`python/swl/transpile/cwl/emit.py`)

The key design question is how to represent a batch workflow in:
1. compiled SWL JSON
2. packed CWL

### Decision summary

Distinguish clearly between task, simple workflow, and batch workflow.

- **task**: shell-script unit, `rec -> rec`
- **simple workflow**: ordinary workflow, `rec -> rec`
- **batch workflow**: workflow whose root parameter is `[rec]`, and whose result is `rec` or `[rec]`

Representation decisions:
- `map` is the source of batch lifting
- mapped callees are still ordinary tasks or ordinary workflows; they are not reclassified as batch callees
- mapping over batch workflows is forbidden for now
- compiled JSON must preserve mapped execution symbolically so that runtime batch size is not required at compile time
- the compiler artifact must contain enough information for the CWL transpiler to lower mapped execution to runnable CWL scatter
- array projection such as `calls.bcf` must be represented explicitly and type-checked against the homogeneous mapped element output schema

### Why this representation

This should follow the actual language structure more closely.

In the example:
- `call_variant.wf` is a simple workflow
- `merge.sh` is a task
- `panel.swl` is the batch workflow
- inside `panel.swl`, the expression `map call_variant xs` means that the imported callee is executed in mapped mode over the batch input

So in compiled form we need to preserve these distinctions:
- a normal task call producing `rec`
- a normal workflow expansion producing `rec`
- a mapped task or mapped workflow step producing `[rec]`
- a field projection from `rec`
- a field projection from `[rec]`

The important new concept is mapped execution. The compiled representation must preserve it explicitly rather than fully expanding it away.

---

## 1. Scope the final target and early implementation boundaries

Final target semantics:
- `map` accepts any function whose effective signature is `rec -> rec`
- `map f` has type `[rec] -> [rec]`
- `.` works on `[rec]` and returns `[t]`
- tasks may consume array inputs if the task annotation type is an array type like `[file]`
- batch workflow root type is `[rec] -> rec` or `[rec] -> [rec]`
- mapped results are homogeneous
- missing projected fields are compile-time errors
- batch workflows can be transpiled to runnable packed CWL using native CWL scatter

Implementation boundaries for the first delivered version may still reject some forms if necessary, but the design target should remain the above.
In particular, the implementation should be structured so that it can support all `rec -> rec` mapped callees, not just imported workflows.

Still deferred for later:
- mapping over batch workflows
- nested arrays beyond the direct result of `map`
- multi-array zip/map semantics beyond what is needed for CWL scatter lowering
- arbitrary synthesis of complex CWL `valueFrom` / `linkMerge` expressions
- runtime/input surface syntax for external users to supply `[rec]` batches

---

## 2. Extend the type system for batch values

Files:
- `python/swl/semantic/task/type.py`
- workflow semantic checker code under `python/swl/semantic/wf/`

Plan:

1. Add workflow-level value categories for:
   - `record`
   - `[record]`
   - scalar arrays like `[file]`, `[str]`, `[int]`, `[float]`
   - function `rec -> rec`
   - function `[rec] -> rec`
   - function `[rec] -> [rec]`

2. Keep task annotation array types as-is (`[file]`, `[str]`, `[int]`, `[float]`), but add workflow reasoning for collection-valued expressions and batch workflow roots.

3. Add type rules:
   - if `f : rec -> rec`, then `map f : [rec] -> [rec]`
   - if `xs : [rec]` and `field : t` is present in the element record type, then `xs.field : [t]`
   - if `field` is absent from the element record type, projection is a compile-time error
   - if `f` is batch-typed, then `map f` is a compile-time error for now

4. Ensure workflow well-formedness distinguishes:
   - simple workflow root signature: `rec -> rec`
   - batch workflow root signature: `[rec] -> rec` or `[rec] -> [rec]`

Recommendation:
- do not try to force this into `TaskSignature` alone
- add a workflow type layer for workflow expressions, including element schemas for `[rec]`
- keep mapped record arrays homogeneous by construction

---

## 3. Add explicit IR nodes for batch semantics

Files:
- `python/swl/ir/node.py`
- `python/swl/ir/lower.py`

Plan:

1. Add IR nodes such as:
   - `Map(function, arg)` or equivalent builtin-lowered form
   - `ArrayField(record_array, name)` for projection from `[rec]` to `[t]`

2. Recognize builtin `map` in lowering.
   - process `map` through the same builtin-dispatch path as `import`
   - lower `map f xs` into dedicated IR rather than leaving it as ordinary generic apply

3. Extend field access lowering/typing so `calls.bcf` can refer either to:
   - `ir.Field(record, "bcf")`
   - `ir.ArrayField(record_array, "bcf")`

4. Carry enough type/schema information to distinguish:
   - ordinary record field projection
   - array-of-record field projection
   - compile-time invalid projection on missing element fields

Recommendation:
- prefer explicit IR nodes over silently overloading existing `Field` semantics
- make builtin recognition explicit early so that later forcing and transpilation can rely on structured IR

---

## 4. Extend forced DAG values and compiled JSON

Files:
- `python/swl/ir/dag.py`
- `python/swl/ir/force.py`

This is the main representation decision.

### 4.1 New forced value kinds

Add forced/runtime-planning value kinds in `dag.py` for mapped execution and array projection.
For example:
- an explicit mapped call/step value
- an explicit array-field projection value
- array input metadata sufficient for batch typing and later CWL lowering

Suggested in-memory shape:

```python
@dataclass
class MappedCall:
    id: str
    path: str
    callee_kind: str            # 'task' or 'workflow'
    source: object              # array-valued driver
    outputs: List[str]          # element-record output fields
    input_schema: Optional[dict] = None
    output_schema: Optional[dict] = None
    deps: List[str] = field(default_factory=list)
```

and:

```python
@dataclass(frozen=True)
class ArrayField:
    source: object
    name: str
```

Whether the top-level compiled JSON keeps `tasks[]` or later generalizes to `steps[]`, the artifact must preserve mapped execution symbolically and encode enough callee metadata for CWL lowering.

### 4.2 JSON representation requirements

The compiled artifact must support all of the following concepts:
- ordinary concrete task call
- mapped task call
- mapped workflow call
- projection from mapped `[rec]` result to array output like `[file]`
- top-level workflow outputs that may be `rec` or `[rec]`

Two schema directions are possible:

#### Option A: generalize `tasks[]` to `steps[]`

This gives a uniform representation where each executable unit is either a task or workflow step, and mapped execution is an explicit property.

Example:

```json
{
  "inputs": {
    "xs": {"type": "[rec]"},
    "outbase": {"type": "str"}
  },
  "steps": [
    {
      "id": "call_variant",
      "type": "workflow",
      "path": "call_variant.wf",
      "map": {
        "source": {"source": "input", "name": "xs"}
      },
      "deps": [],
      "outputs": {
        "bcf": {"type": "file"}
      }
    },
    {
      "id": "merge",
      "type": "task",
      "path": "merge.sh",
      "deps": ["call_variant"],
      "bindings": {
        "bcf": {
          "source": "array_field",
          "step": "call_variant",
          "output": "bcf"
        }
      }
    }
  ],
  "outputs": {
    "bcf": {"source": "step", "step": "merge", "output": "bcf"}
  }
}
```

#### Option B: preserve `tasks[]` and add explicit mapped/task-or-workflow variants

This reduces schema churn for the current codebase, but the representation must still encode mapped workflow or mapped task execution explicitly.
A task-only schema is not sufficient for the final target, because runnable CWL transpilation needs to know whether a mapped callee lowers to a `CommandLineTool` or a subworkflow.

### 4.3 Practical recommendation

For the final target, a uniform executable collection is cleaner because batch lowering needs to represent mapped tasks and mapped workflows explicitly.
However, implementation may stage toward that if needed.

Required properties regardless of exact field names:
- callee kind must be explicit (`task` vs `workflow`)
- mapped execution must be explicit
- array-field projection must be explicit
- references must preserve dependencies and output schemas
- the artifact must carry enough information for the CWL emitter to produce runnable scatter steps

### 4.4 Why not pre-expand the batch into many ordinary calls?

Do not compile batch execution by cloning the mapped workflow into many ordinary steps in the JSON.

Reasons:
- batch cardinality is runtime data-dependent
- compiled JSON should stay symbolic and self-contained
- CWL scatter is already symbolic and array-driven
- pre-expansion would require concrete batch size at compile time

So the compiled artifact should represent intent: “this call is scattered over this array input”.

---

## 5. Forcing semantics for `map`

File:
- `python/swl/ir/force.py`

Plan:

1. Teach `force_value` / `_apply` to evaluate `ir.Map`.

2. When forcing `map f xs`:
   - require `f` to have effective signature `rec -> rec`
   - reject `f` if it is batch-typed
   - require `xs` to be an array-of-record value with known homogeneous element schema
   - produce a mapped symbolic value, not a concrete list of calls

3. When forcing field projection from a mapped result:
   - `ArrayField(MappedCall, "bcf")` should represent `[file]`
   - if `bcf` is not in the mapped element output schema, reject at compile time

4. Extend dependency walking:
   - downstream consumers of `calls.bcf` should depend on the mapped callee step

5. Extend output normalization and serialization/deserialization for new sources:
   - array-field projection
   - mapped step references
   - top-level `[rec]` outputs where supported

6. Ensure forcing can preserve enough callee metadata for CWL lowering:
   - task vs workflow
   - input signature
   - output signature
   - mapped source wiring

Recommendation:
- forcing should create a symbolic mapped step/call representation rather than expanding the batch
- this keeps runtime cardinality symbolic and aligns with later CWL scatter lowering

---

## 6. CWL lowering target

File:
- `python/swl/transpile/cwl/emit.py`

Batch workflows should transpile to runnable packed CWL using native CWL `Workflow` scatter.

### 6.1 Main lowering rule

For:

```swl
calls = map call_variant xs
merge { bcf: calls.bcf }
```

the transpiler should emit a scattered workflow step for `calls`, and then a downstream reduction step consuming the scattered array output.

Mapped task callees should lower to scattered `CommandLineTool` steps.
Mapped workflow callees should lower to scattered subworkflow steps.

### 6.2 Important CWL decision

Do not represent batch structure in CWL by inventing a custom runtime convention.
Use native CWL scatter.

Reason:
- scatter is the direct analog of SWL `map`
- it preserves dependency semantics
- it gives portable execution behavior in CWL runners
- it is the path to runnable CWL

### 6.3 CWL type rules

Extend `_cwl_type` to support:
- `[file]` -> `{ "type": "array", "items": "File" }`
- `[str]` -> `{ "type": "array", "items": "string" }`
- `[int]` -> `{ "type": "array", "items": "int" }`
- `[float]` -> `{ "type": "array", "items": "float" }`

The transpiler must also be able to lower mapped `[rec]` flow into CWL-compatible scattered ports.

### 6.4 Resolve the `[rec]` to CWL mismatch

CWL supports records and arrays, but the current compiler does not yet preserve enough record-array schema information in its DAG artifact.
So the final target requires an explicit normalization path for CWL lowering.

Recommended final approach:

1. The canonical compiler artifact may preserve the logical `[rec]` semantics.
2. Before emitting CWL scatter, the transpiler must normalize the mapped `[rec]` source into concrete scattered ports using the mapped callee input schema.

Example conceptual rewrite:
- logical batch input: `xs : [rec{fastq1:file, fastq2:file, ref:file, outbase:str}]`
- normalized scattered ports:
  - `fastq1 : File[]`
  - `fastq2 : File[]`
  - `ref : File[]` or scalar `File` if bound outside the mapped item
  - `outbase : string[]`

Then the scattered CWL step wires each required input from the appropriate source and uses `scatterMethod: dotproduct`.

This means the transpiler needs access to:
- the mapped callee input signature
- which inputs are driven per-item from the batch element record
- which inputs are pre-bound constants or outer workflow values

### 6.5 Runnable CWL requirement

The final target is not just symbolic JSON.
The design must ensure that a batch workflow can be transpiled into runnable packed CWL.
That implies:
- mapped tasks lower to runnable scattered tool steps
- mapped workflows lower to runnable scattered subworkflow steps
- array field projection such as `calls.bcf` becomes an array-valued CWL source
- downstream reduction tasks can consume that array-valued source
- `[rec] -> [rec]` workflows can either lower to scattered step outputs exposed as workflow outputs, or be rejected only where the target runner format truly cannot express them yet

If some subcase remains unsupported during phased implementation, it should be rejected explicitly, but the architecture should still target runnable CWL for the supported batch cases.

---

## 7. Concrete phased implementation steps

### Phase 1: workflow typing and IR

1. Add workflow-level typing for:
   - `rec`
   - `[rec]`
   - scalar arrays
   - function categories `rec -> rec`, `[rec] -> rec`, `[rec] -> [rec]`
2. Add compile-time checking for:
   - `map` only on effective `rec -> rec`
   - rejection of `map` on batch workflows
   - `.` on `[rec]` with missing fields rejected at compile time
3. Add IR support for:
   - `ir.Map`
   - `ir.ArrayField`
4. Route builtin `map` through the same builtin-recognition path used for `import`

Deliverable:
- semantic checking and lowering recognize batch constructs and reject invalid ones correctly.

### Phase 2: forced symbolic batch representation

1. Add forced DAG/runtime value classes for:
   - mapped call/step
   - array-field projection
   - array-aware input/output metadata
2. Extend serialization/deserialization with explicit mapped and array-field sources.
3. Preserve enough metadata for downstream CWL normalization:
   - callee kind
   - path
   - input schema
   - output schema
   - mapped source
4. Add tests for forcing batch workflows into stable symbolic JSON.

Deliverable:
- `python/swl/eval_force.py` on batch workflows prints stable symbolic JSON that preserves mapped execution.

### Phase 3: CWL normalization and scatter emission

1. Extend emitter validation to accept mapped steps and array-field bindings.
2. Add array type lowering in `_cwl_type`.
3. Implement normalization from logical `[rec]` mapping to concrete scattered CWL ports.
4. Emit runnable scattered CWL steps:
   - scattered tool step for mapped task
   - scattered subworkflow step for mapped workflow
   - `scatterMethod: dotproduct`
5. Wire downstream reduction steps from scattered array outputs.
6. Add explicit rejection for any still-unsupported corner cases.

Deliverable:
- supported batch workflows transpile to runnable packed CWL.

### Phase 4: top-level batch outputs and polish

1. Ensure `[rec] -> [rec]` workflows can be emitted as runnable CWL where representable.
2. Add end-to-end golden tests from `.swl` -> DAG JSON -> packed CWL.
3. Update `spec.md` later with batch syntax/semantics, including `map`.
4. Document remaining deferred topics such as external batch input surface syntax.

Deliverable:
- final target behavior documented and covered by tests.

---

## 8. Test plan

Add tests under:
- `python/swl/ir/test_force.py`
- `python/swl/ir/test_force_codec.py`
- `python/swl/transpile/cwl/test_emit.py`
- workflow semantic tests under `python/swl/semantic/wf/`

Minimum cases:

1. `map` over imported simple workflow produces one mapped symbolic callee in compiled JSON.
2. `map` over imported task also produces one mapped symbolic callee in compiled JSON.
3. `map` over lambda or partial function with effective signature `rec -> rec` type-checks and lowers correctly.
4. `map` over batch workflow is rejected.
5. `calls.bcf` serializes as explicit array-field projection.
6. missing projected field on `[rec]` is rejected at compile time.
7. downstream task consuming `[file]` from `calls.bcf` gets correct dependency.
8. batch workflow returning `rec` transpiles to runnable CWL with scatter + reduction.
9. batch workflow returning `[rec]` transpiles where supported, or is rejected explicitly and locally if not yet implemented.
10. unsupported nested batch cases are rejected explicitly.

---

## 9. Recommended final decisions

### Language and typing

- `map` is a builtin handled in the same builtin-processing path as `import`
- `map` may be applied to any function with effective signature `rec -> rec`
- `map` on batch workflows is forbidden for now
- `[rec]` produced by `map` is homogeneous
- `.` on `[rec]` returns `[t]`
- missing projected fields are compile-time errors
- batch workflow root signatures are `[rec] -> rec` and `[rec] -> [rec]`

### Compiled artifact

- preserve mapped execution symbolically
- make callee kind explicit (`task` vs `workflow`)
- represent array-field projection explicitly
- preserve enough schema/signature information for target-specific lowering
- do not pre-expand batch cardinality into concrete repeated calls

### CWL

- transpile batch workflows using native CWL scatter
- normalize logical `[rec]` mapping into concrete scattered ports before emission
- support mapped tasks and mapped workflows
- ensure the supported batch subset transpiles to runnable packed CWL

---

## 10. Immediate next coding tasks

1. implement workflow-level batch typing and batch workflow root classification
2. implement builtin recognition for `map` alongside `import`
3. add `ir.Map` and `ir.ArrayField`
4. add forced mapped-call and array-field representations plus codec support
5. preserve enough callee schema/signature metadata for CWL lowering
6. implement CWL normalization from logical `[rec]` mapping to scattered ports
7. emit runnable scattered CWL for mapped task/workflow calls
8. add focused force/codec/CWL tests for both `[rec] -> rec` and `[rec] -> [rec]` workflows

That should make the final target explicit while keeping deferred items clearly separated.