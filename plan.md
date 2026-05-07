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

This gives clear boundaries between:
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
- runs task syntax diagnostics on all `tests/*.sh`
- runs task semantic diagnostics on all `tests/*.sh`
- runs workflow syntax diagnostics/evaluation on all `tests/*.swl`
- runs workflow semantic diagnostics on all `tests/*.swl`

#### Unit tests added
- `python/swl/syntax/task/test_parser.py`
- `python/swl/syntax/task/test_interpolation.py`
- `python/swl/syntax/task/test_bash.py`

These cover the basic parsing surface and current shell examples.

## Phase 3: Semantic layer - FIRST PASS COMPLETE

### 3.1 Task semantic typing
Implemented in:
- `python/swl/semantic/task/type.py`

Current behavior:
- converts parsed task annotations into `TaskSignature`
- preserves defaults as interpolation syntax objects
- checks duplicate names and basic chain compatibility

### 3.2 Workflow semantic checking
Implemented in:
- `python/swl/semantic/wf/check.py`

Current behavior:
- resolves imported `.sh` tasks
- resolves imported `.swl` workflows
- builds signatures for both tasks and workflows
- checks explicit `chain` compatibility
- performs conservative workflow input inference
- performs approximate workflow output/signature inference
- detects circular workflow imports

### 3.3 Semantic diagnostics
Implemented in:
- `python/swl/eval_task_semantic.py`
- `python/swl/eval_wf_semantic.py`

Current behavior:
- task semantic diagnostics print semantic signatures
- workflow semantic diagnostics print:
  - imports
  - import kinds
  - chain errors
  - inferred inputs
  - inferred workflow signature

### 3.4 Current semantic model
The current semantic checker already uses a very small symbolic value model:
- open records
- closed records
- task/workflow results
- unknown values

This has been useful for inference, but it is still only an approximation layer.

## Phase 4: Settle evaluation model before IR - DONE

The core semantic decision for IR is now:

### 4.1 Lazy evaluation
The language uses lazy evaluation.

That means:
- expressions do not eagerly execute tasks/workflows when syntactically applied
- evaluation produces values or thunks/closures in the semantic/runtime model
- actual execution should only be forced when building or realizing the execution graph

### 4.2 Functions everywhere
Imported tasks and imported workflows are both functions.

That means:
- `import "align.sh"` yields a function value
- `import "subworkflow.swl"` yields a function value
- user-defined workflow lambdas are also function values

These should all be represented uniformly in IR.

### 4.3 Partial application semantics
Partial application returns a function.

This is now the intended meaning of expressions like:
- `align_hg38 = align { ref: ..., ref_fai: ... }`

This should **not** be modeled as immediate task execution.
Instead:
- applying a function to some arguments returns a new function when required inputs remain unsatisfied
- applying a function with enough inputs returns a lazy computation value
- forcing that lazy computation later contributes nodes to the execution IR

### 4.4 Consequence for current checker
The current workflow checker does not yet fully implement this lazy partial-application model.
It still approximates such expressions using immediate-demand reasoning.
That is acceptable temporarily, but IR work should replace that approximation with a first-class lazy function model.

## Phase 5: IR design - NEXT

The next major phase is to build a proper intermediate representation that matches the lazy semantics above.

### 5.1 IR goals
The IR should:
- represent lazy evaluation, not eager execution
- represent tasks, imported workflows, and lambdas uniformly as functions
- support partial application as a first-class operation
- preserve enough structure for type/input/output checking
- lower naturally into an execution DAG when computation is forced

### 5.2 Proposed runtime/IR value kinds
Suggested value categories:

#### Primitive / literal values
- strings
- ints
- floats
- possibly booleans later

#### Record values
- finite mappings from field name to value
- updates/merges produce derived record values
- field access projects from record values lazily

#### Function values
A uniform callable value with variants:
- imported task function
- imported workflow function
- lambda function
- partially applied function / closure

Each function value should carry:
- parameter/interface information
- captured environment if needed
- provenance (task/workflow/lambda)

#### Computation values
A lazy application result representing:
- function applied to argument(s)
- not yet executed
- available for further composition, projection, or forcing

This is the key value kind missing today.

### 5.3 Proposed IR nodes
At minimum, introduce explicit IR nodes roughly like:
- `IRLiteral`
- `IRRecord`
- `IRField`
- `IRUpdate`
- `IRImportTask`
- `IRImportWorkflow`
- `IRLambda`
- `IRClosure`
- `IRApply`
- `IRBind`
- `IRBlock`
- `IRChain`

A second stage may then lower these into execution-oriented nodes such as:
- `ExecTask`
- `ExecWorkflowCall`
- `ExecValue`
- `ExecProjection`
- `ExecMerge`

### 5.4 Two-level IR recommendation
Recommended split:

#### Semantic IR
A lazy functional IR that mirrors source semantics.
Responsibilities:
- name resolution
- lexical binding
- imports as function values
- partial application
- lazy application
- record construction/update/projection

#### Execution IR / DAG
A lower-level forced computation graph.
Responsibilities:
- concrete task/workflow call nodes
- data dependencies
- named outputs
- runtime parameter propagation
- eventual scheduling/execution

This split should reduce semantic confusion.
Do not make the first IR immediately execution-specific.

### 5.5 Suggested function model
Define a single callable abstraction, e.g. conceptually:
- `FunctionValue(kind, signature, env, body_or_target, bound_args)`

Where:
- `kind` is one of `task`, `workflow`, `lambda`, `partial`
- `signature` describes required/available inputs and outputs
- `env` captures lexical environment if needed
- `body_or_target` points to task metadata, workflow IR, or lambda body
- `bound_args` stores already supplied arguments

Then application behaves as:
1. merge supplied argument with already bound arguments
2. if required inputs remain missing, return another `FunctionValue` (partial)
3. if sufficiently saturated, return a lazy `ComputationValue`

This directly matches the chosen semantics.

### 5.6 Argument model
Callable argument semantics are now clarified.

Rules:
- workflow/task functions conceptually consume records
- applying a task to a scalar argument is legal
- a scalar argument is lifted to a record with a single field named as the first input of the task
- partial application creates a closure that encloses the provided inputs

Consequences for IR:
- semantic lowering must support scalar-to-record lifting for task application
- the function application path must be able to produce closures
- closure values must remember already-bound inputs for later saturation

### 5.7 Import lowering model
Recommended import behavior:
- `import "task.sh"` lowers to `IRImportTask(path, signature)`
- `import "workflow.swl"` lowers to `IRImportWorkflow(path, signature, referenced_ir?)`

Both evaluate to function values.
They should not be executed by import itself.

### 5.8 Chain lowering model
For `a | b`:
- first resolve both sides as function values
- semantically this is function composition over task/workflow record interfaces
- chain checking remains signature-based
- chain lowering should produce either:
  - a composed function value, or
  - syntactic sugar over lambda/application in semantic IR

Recommendation:
- treat chain as sugar for composition at the semantic IR level
- lower it before execution DAG construction

### 5.9 Interpolation placement in IR
Interpolation timing is now clarified.

Keep task interpolation defaults as structured values until later phases.
Their meanings are:
- `Word`: compile-time verbatim substitution, with compile-time syntax checks where possible
- `Var`: runtime substitution from the runtime value environment
- `Expr`: runtime evaluation followed by substitution of the resulting value

Resolution should happen when:
- a concrete task call node exists
- a concrete environment of bound inputs/run params exists
- pre-runtime bash validation is being performed

This implies two validation opportunities:
- pre-runtime validation after compile-time interpolation pieces are known
- runtime validation after runtime substitutions are known

### 5.10 Output/signature policy
Workflow signature semantics are now clarified.

Rules:
- a workflow must evaluate to a function
- the workflow output is the output of that function
- if the workflow final value is an explicit lambda, the outputs are determined by the final expression in the lambda body
- if the workflow final value is a named task, the workflow output is the task output
- if the workflow final value is a chain, semantic analysis should take the union of the output variables of each task/workflow in the chain from left to right

Task output policy:
- task output params must have defaults
- output defaults may include glob patterns such as `*`

Typing policy:
- preserve output types whenever recoverable from task/workflow signatures
- record outputs synthesized by workflow analysis may still have unknown type where necessary

## Phase 6: Concrete implementation steps for IR

### 6.1 Add semantic IR module
Suggested new package:
- `python/swl/ir/`

Suggested files:
- `python/swl/ir/node.py`
- `python/swl/ir/lower.py`
- `python/swl/ir/eval.py` or `python/swl/ir/force.py`

### 6.2 First IR milestone: lazy semantic lowering
Implement lowering from workflow AST + imports into semantic IR.

Scope:
- imports
- lambdas
- records
- field access
- update
- application
- chain
- block/bind

No execution yet, only semantic lowering.

### 6.3 Second IR milestone: function values and partial application
Implement function-value objects and application rules:
- imported task as function
- imported workflow as function
- lambda as function
- partial application returns function
- saturated application returns lazy computation value

### 6.4 Third IR milestone: force/lower to execution DAG
Only after semantic IR is stable:
- force selected computation roots
- build task/workflow call graph
- propagate concrete inputs
- resolve interpolation for task execution nodes

### 6.5 Diagnostics for IR
Add debugging entrypoints such as:
- `python -m swl.eval_ir file.swl`
- print semantic IR
- optionally print forced execution DAG

## Phase 7: What still remains non-blocking

These are important, but they should not block starting IR:
- implementing the two-stage bash validation promised by the clarified interpolation model
- perfect workflow type inference
- precise provenance tracking for every field
- preserving exact output types for every inferred workflow record

## Recommended immediate next step

Implement the first semantic IR with the lazy function model:
- imports produce function values
- lambdas produce function values
- partial application returns closure/function values
- saturated application returns lazy computation values
- chain lowers to function composition
- task application supports scalar-to-record lifting using the first declared input name
- task outputs require defaults and may include glob patterns

That is the cleanest bridge from the current checker to a real execution model.
