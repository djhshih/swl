# Implementation Plan

## Phase 1: Workflow parser - COMPLETE

The workflow parser is structured as:
- `python/swl/syntax/wf/lexer.py`
- `python/swl/syntax/wf/node.py`
- `python/swl/syntax/wf/parser.py`

This remains the model for workflow syntax.

## Phase 2: Task syntax refactor - MOSTLY COMPLETE

### Current implemented structure

Task syntax is now organized as:
- `python/swl/syntax/task/node.py`
- `python/swl/syntax/task/parser.py`
- `python/swl/syntax/task/interpolation.py`
- `python/swl/syntax/task/bash.py`
- `python/swl/eval_task.py`

Semantic task typing is now separated as:
- `python/swl/semantic/task/type.py`

This is a major improvement over the old state because we now have clear boundaries between:
- task annotation parsing
- interpolation parsing
- optional bash-body analysis
- semantic typing/checking

### What is complete

#### Task annotation parser
Implemented in:
- `python/swl/syntax/task/parser.py`

Current behavior:
- splits task file into annotation region and raw body
- parses task doc line
- parses sections: `in`, `out`, `run`
- parses parameter names, optional type, default, and description
- supports multiline description continuations beginning with `|`
- parses defaults through `task/interpolation.py`
- returns `Task(annotation, body)`

#### Task syntax/data nodes
Implemented in:
- `python/swl/syntax/task/node.py`

Current nodes:
- `Task`
- `Annotation`
- `Section`
- `Param`
- `SectionType`

#### Shared interpolation parser
Implemented in:
- `python/swl/syntax/task/interpolation.py`

Current behavior:
- parses plain literals
- parses `$x`
- parses `${x}`
- parses mixed forms like `${outbase}.bam`
- parses conservative expression forms like `${memory / cpu}`

This parser is shared across:
- annotation defaults now
- bash-body analysis later or optionally now

#### Optional bash-body analyzer
Implemented in:
- `python/swl/syntax/task/bash.py`

Current behavior:
- keeps bash analysis separate from task annotation parsing
- parses simple assignments conservatively
- parses command lines conservatively
- extracts interpolation-bearing shell words using the shared interpolation parser

Important design choice:
- `Task` still stores raw body text only
- `bash.py` remains a separate analysis step
- this is intentional and keeps the task parser simpler

#### Diagnostic parser entrypoint
Implemented in:
- `python/swl/eval_task.py`

This mirrors `python/swl/eval.py` and prints:
- annotation doc
- sections
- params
- defaults
- body

#### Test integration
Implemented in:
- `test.sh`

Current behavior:
- runs unit tests
- runs task parser diagnostics on all `tests/*.sh`
- runs workflow diagnostics/evaluation on all `tests/*.swl`

#### Unit tests added
- `python/swl/syntax/task/test_parser.py`
- `python/swl/syntax/task/test_interpolation.py`
- `python/swl/syntax/task/test_bash.py`

These cover the basic parsing surface and current shell examples.

## Phase 3: Syntax cleanup before semantics - RECOMMENDED

Before working deeply on semantics and interpretation, there are a few worthwhile cleanup steps.

These are not architectural rewrites, just polishing the syntax layer so the semantic layer has stable inputs.

### 3.1 Normalize task annotation shape for downstream use

Right now `Annotation` stores `sections`, and downstream code will probably repeatedly search them.

Possible improvement:
- keep `sections` for preserving source structure
- also add convenient accessors or normalized views:
  - `annotation.inputs`
  - `annotation.outputs`
  - `annotation.run`

This is optional, but likely useful before semantic work.

### 3.2 Decide how strict the task parser should be

Open questions to settle before semantics:
- should duplicate names in one section be a parse error or semantic error?
- should duplicate names across `in/out/run` be allowed syntactically?
- should empty sections be allowed?
- should task files require exactly one doc line?
- should body-less tasks be allowed?

Recommendation:
- keep parser responsible only for surface syntax
- move most name-validation and consistency checks into semantics
- but explicitly document these boundaries now

### 3.3 Add a few missing parser tests

Useful remaining tests:
- optional types like `file?`
- array types like `[file]`
- default values with quotes
- duplicate section headers
- empty sections
- malformed default text in annotations
- `$x` form in interpolation, not just `${x}`
- interpolation with multiple literals and vars combined

### 3.4 Decide whether `bash.py` is actually needed

Current status:
- it exists
- it is separate
- it is intentionally conservative

Recommendation:
- keep it, but treat it as provisional
- do not let it drive semantic design yet
- if semantics can proceed from annotation metadata alone, do that first
- only deepen bash analysis when interpretation actually needs it

This matches the current preference:
- raw body remains in `Task`
- bash analysis remains optional and separate

## Phase 4: Semantic layer design - NEXT

The next major phase should be semantics, built on top of the stabilized syntax layer.

### 4.1 Task signature construction

Use parsed task annotations to construct semantic task signatures.

Likely input:
- `syntax.task.node.Task`

Likely output:
- `semantic.task.type.TaskSignature`

Responsibilities:
- convert annotation parameter types into `TypeKind`
- convert `in/out/run` sections into normalized maps
- preserve defaults for later evaluation/substitution

Suggested API:
- `signature_from_task(task) -> TaskSignature`

This is the bridge between syntax and semantics.

### 4.2 Validation of task signatures

Semantic validation should likely include:
- duplicate parameter names within a section
- duplicate parameter names across sections if disallowed
- output params should probably have types
- required vs optional input interpretation
- consistency of runtime params if any constraints exist

These are better handled in semantics than parsing.

### 4.3 Workflow semantic checking

Then build workflow semantics on top of:
- parsed workflow AST
- imported task signatures

This phase should handle:
- import resolution
- type compatibility in chains
- compatibility in record update / merge where relevant
- workflow input inference
- workflow output inference if needed
- circular import checks
- DAG construction and cycle detection

### 4.4 Interpretation / evaluation

Only after task signatures and workflow semantics are stable should we move into interpretation.

Likely steps:
- `import "task.sh"` returns a callable task value/signature
- `import "workflow.swl"` returns a callable workflow value
- evaluate workflow AST into an executable DAG or intermediate form
- perform partial application and chaining semantics
- resolve task defaults and input propagation

## Phase 5: What still may be needed before interpretation

Strictly speaking, we can start semantics now.

But a few small things may still be helpful first.

### 5.1 Add syntax-to-semantic bridge functions

This is probably the single most useful next step before broader interpretation.

Examples:
- parse task file -> `Task`
- build task signature -> `TaskSignature`
- parse workflow file -> workflow AST

This avoids mixing parsing with semantic logic.

### 5.2 Decide how interpolation should be represented semantically

Question:
- should defaults remain syntax nodes (`Word`, `Var`, `Expr`) until execution?
- or should some be normalized earlier?

Recommendation:
- keep interpolation as syntax objects for now
- only resolve them during later evaluation when a value environment exists

That keeps semantics simpler.

### 5.3 Clarify how much shell semantics we need

Important question before interpretation:
- do workflow/task semantics depend only on annotations?
- or do we need to understand the bash body to determine outputs or dependencies?

Given current examples, annotations already declare inputs/outputs/run.

Recommendation:
- treat annotation metadata as authoritative for semantics
- treat bash body as opaque execution text for now
- only use `bash.py` later for optional validation or interpolation support if needed

This means we probably do not need deeper shell parsing before starting semantics.

## Recommended next steps

### Immediate next step
Implement a syntax-to-semantic bridge for tasks:
- `Task` -> `TaskSignature`

That will make the task parser immediately useful to the semantic layer.

### Then
Implement semantic import handling and workflow checking:
- resolve imported `.sh` and `.swl` files
- construct signatures for imported tasks
- check type compatibility in workflow chains
- infer workflow inputs

### Then
Implement interpretation:
- evaluate imports into callable values
- evaluate workflow expressions
- build DAG/intermediate execution representation

## Phase 6: Concrete next implementation steps - IN PROGRESS

### 6.1 Build task signature bridge
- add `signature_from_task(task)` in `python/swl/semantic/task/type.py`
- map parsed task annotation params into:
  - `inputs`
  - `outputs`
  - `run`
- parse param type strings with `parse_type()`
- preserve defaults as syntax/interpolation objects for later evaluation
- detect duplicate parameter names at semantic boundary

### 6.2 Add semantic tests for task signatures
- valid task -> `TaskSignature`
- duplicate input/output names fail
- missing required types where disallowed fail
- run params are preserved
- outputs with interpolation defaults are preserved

### 6.3 Add diagnostic entrypoint for semantic task signature
- add something like `python/swl/eval_task_semantic.py`
- print parsed task plus constructed `TaskSignature`
- use this from `test.sh` if helpful

### 6.4 Start workflow semantic import layer - DONE (first pass)
- parse workflow
- find imported task paths
- load imported task files
- build task signatures for imported tasks
- perform first-pass chain checking and workflow input inference

### 6.5 Build workflow semantic IR - IN PROGRESS
- introduce symbolic workflow values for:
  - records
  - task results
  - unknown values
- evaluate workflow bodies conservatively
- support update/record/function workflows better than simple chain-only analysis
- replace ad-hoc input inference with proper demand-driven record-field inference
- build approximate workflow signatures from workflow bodies
- support importing `.swl` workflows semantically, not just `.sh` tasks

### 6.6 Record issues during implementation
- document mismatches between current code, spec, and desired semantics in `issues.md`
- especially note unresolved questions around:
  - interpolation resolution timing
  - whether outputs must always have defaults
  - how much shell analysis is needed for execution

## Summary answer: is there anything else to do before semantics and interpretation?

Not much.

The main parser architecture is now in a good place.

Before semantics, the only really worthwhile remaining syntax-side work is:
- a few more focused parser tests
- possibly small convenience accessors on `Annotation`
- a clear bridge from parsed `Task` to semantic `TaskSignature`

We do **not** need to fully parse bash before starting semantics.

Recommended principle going forward:
- annotations drive semantics
- workflow AST drives composition
- interpolation stays as syntax until evaluation
- raw bash body remains opaque unless execution support later requires deeper analysis
