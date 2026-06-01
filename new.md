# Core decision

Tasks and workflows are record-oriented, but with batch support there are now distinct categories that should be described separately.

- `rec`: `{ str: file|num|str, ... }`
- `tab`: record of arrays; each entry must be an array; each array must be the same length
- **task**: `rec -> rec`
- **simple workflow**: `rec -> rec`
- **batch workflow**: `tab -> rec` or `tab -> tab`

For the first implementation target, `map` creates batch behavior. A mapped callee may be:
- a task (`rec -> rec`)
- a simple workflow (`rec -> rec`)

## 1. Task

A task is the existing shell-script unit described by task annotations.

- signature: `rec -> rec`
- compiled representation: one task definition/call unit
- this does **not** change semantically for the batch feature

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
- semantically, a simple workflow is an ordinary workflow over records
- this does **not** change semantically for the batch feature

## 3. Batch workflow

A batch workflow is a workflow whose parameter is `tab`.
A batch workflow may return either:
- `rec`
- `tab`

Example reduction case (`tab -> rec`):

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

Example passthrough case (`tab -> tab`):

```panel_only.swl
call_variant = import "call_variant.wf"

\xs ->
   map call_variant xs
```

Remarks
- `map` lifts a function from `rec -> rec` to `tab -> tab`
- `panel.swl` is a batch workflow with signature `tab -> rec`
- `panel_only.swl` is a batch workflow with signature `tab -> tab`
- `calls` is a value of type `tab`
- get operator `.` extends to `tab` and returns an array: `calls.bcf` returns `[file]`
- because `map` applies one statically typed `rec -> rec` function uniformly to all elements, `tab` produced by `map` is homogeneous, not heterogeneous
- unlike a simple workflow, a batch workflow must preserve symbolic mapped execution semantics; runtime batch cardinality is not part of the source-language semantics

## New builtins and operators

### `map`

`map` is a builtin.

`map` takes a function `f1` with signature `rec -> rec` and returns a function
`f2` with signature `tab -> tab`.

Allowed mapped callees in the target design:
- imported task
- imported simple workflow
- lambda of type `rec -> rec`
- partially applied function whose remaining signature is `rec -> rec`

Forbidden for now:
- `map` applied to a batch workflow (`tab -> rec` or `tab -> tab`)
- nested batch mapping where the mapped callee itself consumes `tab`

When `f2` is applied to a table `tab`, `f1` is run on each logical row of the table,
producing output records that are then reassembled into a homogeneous table.

### `map_by`

`map_by` is a builtin.

`map_by` takes a function `f1` with signature `rec -> rec`, a string key naming a table column,
and returns a function with signature `tab -> tab`.

Allowed mapped callees in the target design:
- imported task
- imported simple workflow
- lambda of type `rec -> rec`
- partially applied function whose remaining signature is `rec -> rec`

Forbidden for now:
- `map_by` applied to a batch workflow (`tab -> rec` or `tab -> tab`)
- nested batch mapping where the mapped callee itself consumes `tab`

When `map_by f1 key` is applied to a table `tab`:
- rows are partitioned by equality of the values in the column named by `key`
- each partition is presented to `f1` as a grouped slice with record shape
- that grouped slice is a `rec`, not a `tab`
- fields in that grouped slice may contain arrays collected from the rows in the group
- those arrays do not by themselves make the grouped slice a `tab`
- one output record is produced per group
- output rows therefore correspond to unique key values, not original input rows
- the grouping key must already exist in the input table schema and must be preserved in the output

### `.` on `tab`

If `xs : tab` and `field : [t]` is a top-level table column, then:

```swl
xs.field
```

returns `[t]`.

If `field` is absent from the statically known table schema, this is a compile-time error.

Do not add a separate field-access operator for table columns; this remains ordinary field access with type-directed meaning.

## Typing rules

This section makes the intended typing discipline explicit.

### 1. Base value categories

At workflow level, the important value categories are:
- scalar types: `file`, `str`, `int`, `float`
- scalar array types: `[file]`, `[str]`, `[int]`, `[float]`
- `rec`: a record with named fields
- `tab`: a table whose top-level fields are homogeneous arrays of equal length
- function types:
  - `rec -> rec`
  - `tab -> rec`
  - `tab -> tab`

A `tab` is record-shaped, but it is not typed as an ordinary `rec`.
This distinction is needed because:
- field access has different result types on `rec` vs `tab`
- `map` specifically consumes `tab`
- simple and batch workflows are classified by this distinction

So the design should treat `tab` as a distinct workflow type, not merely as a conventionally shaped `rec`.

### 2. Record and table schemas

If:
- `r : rec{ a:t1, b:t2, ... }`
- `xs : tab{ a:[t1], b:[t2], ... }`

then the schema records the statically known fields or columns.

A table schema is homogeneous by construction:
- each column has one element type
- all top-level columns have the same logical row count

Rows are derived semantically by zipping corresponding array positions across columns, but rows are not a separate surface type.

### 3. Field access

Field access is type-directed.

For records:
- if `r : rec{ ..., f:t, ... }`, then `r.f : t`
- if `f` is absent from the statically known record schema, that is a compile-time error unless the surrounding typing rule intentionally permits open-record inference

For tables:
- if `xs : tab{ ..., f:[t], ... }`, then `xs.f : [t]`
- if `f` is absent from the statically known table schema, that is a compile-time error

So `.` is one operator with two related typing rules:
- record projection: `rec -> field`
- table column projection: `tab -> [element]`

### 4. Function application

Application is shape-sensitive.

If:
- `f : rec -> rec`
- `x : rec`

then:
- `f x : rec`

If:
- `g : tab -> rec`
- `xs : tab`

then:
- `g xs : rec`

If:
- `h : tab -> tab`
- `xs : tab`

then:
- `h xs : tab`

Applying a `rec -> rec` function directly to a `tab` is a type error.
Applying a batch function directly to a `rec` is a type error.

Partial application is allowed in the ordinary way, but the remaining function type must still be respected.
In particular, if a partially applied function has effective type `rec -> rec`, then it may be used as the callee of `map`.

### 5. `map`

The fundamental typing rule is:

- if `f : rec -> rec`, then `map f : tab -> tab`

and therefore:

- if `f : rec -> rec` and `xs : tab`, then `map f xs : tab`

More explicitly, if:
- `f : rec{a:t1, b:t2, ...} -> rec{u1:v1, u2:v2, ...}`

then:
- `map f : tab{a:[t1], b:[t2], ...} -> tab{u1:[v1], u2:[v2], ...}`

This explains why `calls.bcf` has type `[file]` when:
- `calls = map call_variant xs`
- `call_variant : rec -> rec{ bcf:file, ... }`

Forbidden typing cases for now:
- if `f : tab -> rec`, then `map f` is a compile-time error
- if `f : tab -> tab`, then `map f` is a compile-time error

So `map` is not a general higher-kinded collection combinator here; it is specifically the lifting from row-wise `rec -> rec` computation to table-wise `tab -> tab` computation.

### 6. `map_by`

The fundamental typing rule is:

- if `f : rec -> rec`, then `map_by f : str -> tab -> tab`

and therefore:

- if `f : rec -> rec`, `key : str`, and `xs : tab`, then `map_by f key xs : tab`

More explicitly, if:
- `xs : tab{ k:[tk], a:[t1], b:[t2], ... }`
- `f : rec{k:tk, a:[t1], b:[t2], ...} -> rec{k:tk, u1:v1, u2:v2, ...}`

then:
- `map_by f "k" xs : tab{ k:[tk], u1:[v1], u2:[v2], ... }`

with one output row per distinct value of `xs.k`.

Additional typing rules:
- the `key` argument must name an existing top-level input table column
- if the named key column is absent from the statically known table schema, that is a compile-time error
- the output of the mapped function must preserve the grouping key
- if the key is not present in the mapped output record schema, that is a compile-time error

As with `map`, forbidden typing cases for now:
- if `f : tab -> rec`, then `map_by f` is a compile-time error
- if `f : tab -> tab`, then `map_by f` is a compile-time error

### 7. Workflow classification

Workflow classification is type-based.

- if the root workflow value has type `rec -> rec`, it is a simple workflow
- if the root workflow value has type `tab -> rec`, it is a batch workflow
- if the root workflow value has type `tab -> tab`, it is a batch workflow

This classification should be determined by the inferred input/output types of the workflow body, not by parameter names such as `x` or `xs`.

### 7. Record update and table update

Record update remains:
- if `r1 : rec` and `r2 : rec`, then `r1 // r2 : rec`
- overlapping fields must have compatible types
- right side overrides left side on overlapping fields

For tables, the intended semantics are stricter:
- `tab // tab` is only valid when overlapping fields remain column-typed and row-compatible
- `tab // rec` and `rec // tab` conceptually broadcast scalar record fields across all rows of the table

However, those table-update rules are part of the intended semantics and may be implemented in phases.
The important typing point is that table-aware update should not be collapsed into ordinary unconstrained record merge semantics.

### 8. Subtyping / relationship between `rec` and `tab`

A `tab` is structurally record-shaped, but it should not be identified with `rec` in the workflow type system.

The intended relationship is:
- `tab` behaves like a specialized record-shaped collection value
- `tab` has its own typing rules
- `tab` is not interchangeable with `rec` for function application

So the design should not rely on “a table is just a record whose fields happen to be arrays”.
That structural observation is true, but insufficient for typing.

## Concrete design target

This section describes the intended language and artifact design, not a file-by-file implementation plan.

### Decision summary

Distinguish clearly between task, simple workflow, and batch workflow.

- **task**: shell-script unit, `rec -> rec`
- **simple workflow**: ordinary workflow, `rec -> rec`
- **batch workflow**: workflow whose root parameter is `tab`, and whose result is `rec` or `tab`

Representation decisions:
- `map` is the source of batch lifting
- mapped callees are still ordinary tasks or ordinary workflows; they are not reclassified as batch callees
- mapping over batch workflows is forbidden for now
- compiled representations must preserve mapped execution symbolically so that runtime batch size is not required at compile time
- artifacts must contain enough information for downstream lowering to preserve mapped execution correctly
- array projection such as `calls.bcf` must be represented explicitly in semantics and type-checked against the homogeneous mapped element output schema

### Why this representation

This follows the actual language structure more closely.

In the example:
- `call_variant.wf` is a simple workflow
- `merge.sh` is a task
- `panel.swl` is the batch workflow
- inside `panel.swl`, the expression `map call_variant xs` means that the imported callee is executed in mapped mode over the batch input

So the design needs to preserve these distinctions:
- a normal task call producing `rec`
- a normal workflow expansion/call producing `rec`
- a mapped task or mapped workflow step producing `tab`
- a field projection from `rec`
- a field projection from `tab`

The important new concept is mapped execution. The representation must preserve it explicitly rather than fully expanding it away.

## Final target semantics

- `map` accepts any function whose effective signature is `rec -> rec`
- `map f` has type `tab -> tab`
- tasks may consume array inputs if the task annotation type is an array type like `[file]`
- batch workflow root type is `tab -> rec` or `tab -> tab`
- mapped results are homogeneous
- missing projected fields are compile-time errors
- batch workflows can be lowered to runnable targets using native scatter/gather-style execution where the backend supports it

Still deferred for later:
- mapping over batch workflows
- nested arrays beyond the direct result of `map`
- multi-array zip/map semantics beyond what is needed for the basic `map` model
- arbitrary synthesis of complex backend-specific expressions
- external surface syntax for supplying batch/table inputs

## Recommended final decisions

### Language and typing

- `map` is a builtin
- `map` may be applied to any function with effective signature `rec -> rec`
- `map` on batch workflows is forbidden for now
- `tab` produced by `map` is homogeneous
- `.` on `tab` returns `[t]` for table columns
- missing projected fields are compile-time errors
- batch workflow root signatures are `tab -> rec` and `tab -> tab`

### Representation

- preserve mapped execution symbolically
- make callee kind explicit (`task` vs `workflow`)
- keep table-column projection explicit in semantics, but do not add a separate field-access operator
- preserve enough schema/signature information for target-specific lowering
- do not pre-expand batch cardinality into concrete repeated calls

### Backends

- lower batch workflows using native scatter/gather mechanisms where the backend supports them
- support mapped tasks and mapped workflows
- ensure the supported batch subset lowers to runnable artifacts

## Immediate next step

Implementation details and file-by-file changes belong in `plan.md`.

`new.md` should remain the design note for:
- the distinction between `rec`, `tab`, task, simple workflow, and batch workflow
- the semantics of `map`
- the meaning of `.` on `tab`
- the requirement that mapped execution remain symbolic
