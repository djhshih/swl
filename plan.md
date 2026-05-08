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
- enforces duplicate-binding scope rules semantically
- permits shadowing from any nested scope into outer scopes
- enforces that a workflow must evaluate to a function

### 3.3 Semantic diagnostics
Implemented in:
- `python/swl/eval_task_semantic.py`
- `python/swl/eval_wf_semantic.py`

Current behavior:
- task semantic diagnostics print semantic signatures
- workflow semantic diagnostics print:
  - imports
  - import kinds
  - chain errors / semantic errors
  - inferred inputs
  - inferred workflow signature

### 3.4 Current semantic model
The current semantic checker already uses a small symbolic value model:
- open records
- closed records
- function values
- closure values
- computation values
- unknown values

This has been useful for inference, but it is still only an approximation layer.
Task semantic scope rules are also now enforced at signature-building time:
- inputs, outputs, and run parameters each have their own scope
- duplicates are rejected within a scope
- input/output name overlap is allowed

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

## Phase 5: IR design - IN PROGRESS

The next major phase is to build a proper intermediate representation that matches the lazy semantics above.

### 5.1 IR goals
The IR should:
- represent lazy evaluation, not eager execution
- represent tasks, imported workflows, and lambdas uniformly as functions
- support partial application as a first-class operation
- preserve enough structure for type/input/output checking
- lower naturally into an execution DAG when computation is forced

### 5.2 Runtime/IR value kinds
Current intended value categories:

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

### 5.3 Current semantic IR nodes
Current semantic IR is implemented under:
- `python/swl/ir/node.py`
- `python/swl/ir/lower.py`
- `python/swl/eval_ir.py`

Current nodes include:
- `Literal`
- `Unknown`
- `Name`
- `Record`
- `Field`
- `Update`
- `Function`
- `Lambda`
- `Closure`
- `Apply`
- `Chain`
- `Bind`
- `Block`

Current lowering decisions:
- import bindings are reduced at compile time into `Function` values
- import-only binds are not emitted as IR `Bind`s
- empty top-level import-only blocks reduce directly to their result
- imported workflow functions carry cached lowered bodies
- imported task/workflow functions are cached by function name in the lowerer
- `eval_ir.py` prints the function cache first, then the lowered tree
- `IRChain`
- `IRName`

A second stage may then lower these into execution-oriented nodes such as:
- `ExecTask`
- `ExecWorkflowCall`
- `ExecValue`
- `ExecProjection`
- `ExecMerge`
- `ExecEdge`

Use a single `IRImport` node rather than separate `ImportTask` and `ImportWorkflow` nodes unless execution forcing proves that separate node classes materially simplify the design. The common semantics at the semantic IR layer are:
- imported thing is a function value
- it has a path
- it has a signature
- it has an import kind (`task` or `workflow`)

So the semantic IR can likely model imports as one node with a `kind` field, and defer task/workflow operational differences to forcing.

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
- `import "task.sh"` lowers to `IRImport(kind='task', path=..., signature=...)`
- `import "workflow.swl"` lowers to `IRImport(kind='workflow', path=..., signature=...)`

Both evaluate to function values.
They should not be executed by import itself.

Why prefer one `IRImport` node over separate `ImportTask` / `ImportWorkflow` nodes?
- at the semantic IR layer they are both just imported callable values
- both participate in application, closure formation, and chaining uniformly
- both need the same core metadata: path, signature, kind
- keeping one node reduces branching in lowering and semantic passes

When separate behavior matters:
- forcing to execution DAG may still branch on `kind`
- task forcing may produce a concrete task-execution node
- workflow forcing may recurse into another lowered workflow body or call boundary

So the recommended split is:
- single import node in semantic IR
- differentiated behavior in `ir/force.py`

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
- `python/swl/ir/force.py`

#### Overall IR shape
Use two IR shapes at two different phases:
- semantic IR should be tree-shaped
- forced execution IR should be graph-shaped / DAG-like

Rationale:
- the source language is naturally tree-shaped: lambda, apply, update, record, field access, block
- tree IR is simpler for lexical scope, lazy semantics, and debugging
- execution naturally wants sharing and explicit dependencies, which is graph-shaped
- therefore: semantic tree first, execution DAG second
Recommended responsibilities:
- define IR node/value classes
- define small enums/tags where useful
- preserve enough metadata for later forcing into an execution DAG
- keep nodes simple and explicit rather than clever

Recommended node families:

##### Primitive/value nodes
- `Literal(value)`
  - for strings, ints, floats
- `Unknown()`
  - only if needed for conservative lowering/debugging

##### Record nodes
- `Record(fields: Dict[str, IRNode])`
  - lazy record construction
- `Field(record: IRNode, name: str)`
  - lazy field projection
- `Update(left: IRNode, right: IRNode)`
  - record merge/update

##### Binding/control nodes
- `Bind(name: str, value: IRNode)`
  - a single lexical binding
- `Block(bindings: List[Bind], result: IRNode)`
  - explicit lexical scope with final result

##### Function nodes
- `Lambda(param: str, body: IRNode)`
  - user-defined workflow lambda
- `ImportTask(name: str, path: str, signature: TaskSignature)`
  - imported task as a function value
- `ImportWorkflow(name: str, path: str, signature: TaskSignature)`
  - imported workflow as a function value
- `Closure(function: IRNode, bound_arg: IRNode)` or `Closure(function: IRNode, bound_fields: Dict[str, IRNode])`
  - partial application result

##### Application/composition nodes
- `Apply(function: IRNode, arg: IRNode)`
  - lazy application, not execution
- `Chain(items: List[IRNode])`
  - explicit composition node initially, even if later desugared

##### Name/reference nodes
- `Name(name: str)`
  - local lexical reference after lowering
  - optional depending on lowering style

##### Optional annotation/signature nodes or metadata
- all callable nodes should carry or be accompanied by signature metadata where available
- record-like nodes may optionally carry inferred field/type metadata later

Important design guidance for `node.py`:
- keep it as a pure representation layer
- no filesystem access
- no parsing
- no import resolution
- no execution decisions
- ideally dataclass-style immutable/simple nodes

#### `python/swl/ir/lower.py`
This file should translate workflow syntax + semantic import information into the semantic IR.

Recommended responsibilities:
- load/check workflow through the existing semantic layer as needed
- lower workflow AST nodes into IR nodes
- resolve imported names to `ImportTask` / `ImportWorkflow`
- convert lexical bindings into explicit IR `Block` / `Bind`
- lower chain syntax into either:
  - explicit `Chain` nodes first, or
  - desugared `Apply`/`Lambda` composition later
- attach signature metadata from task/workflow imports
- preserve lazy semantics: lowering should never execute tasks/workflows

Recommended APIs:
- `lower_file(path: str) -> IRNode`
- `lower_tree(tree, imports) -> IRNode`
- `lower_expr(expr, env, imports) -> IRNode`

Suggested lowering rules:

##### Imports
Workflow syntax:
- `x = import "align.sh"`

Lower to:
- bind `x` to `ImportTask(...)` or `ImportWorkflow(...)`

##### Lambda
Workflow syntax:
- `\x -> body`

Lower to:
- `Lambda("x", lower(body))`

##### Block
Workflow syntax block with intermediate bindings and final expr

Lower to:
- `Block([...bindings...], result)`

##### Record
Workflow syntax:
- `{a: x, b: y}`

Lower to:
- `Record({"a": lower(x), "b": lower(y)})`

##### Field access
Workflow syntax:
- `x.foo`

Lower to:
- `Field(lower(x), "foo")`

##### Update
Workflow syntax:
- `a // b`

Lower to:
- `Update(lower(a), lower(b))`

##### Application
Workflow syntax:
- `f x`

Lower to:
- `Apply(lower(f), lower(x))`

Do not decide at lowering time whether this is:
- a full application
- a partial application
- a scalar-lifted application

Those are semantic/runtime questions for later evaluation/forcing passes.
Lowering should preserve structure, not collapse it.

##### Chain
Workflow syntax:
- `a | b | c`

Initial recommendation:
- lower to `Chain([lower(a), lower(b), lower(c)])`

Reason:
- preserves user intent clearly
- lets later semantic/forcing passes decide whether to interpret chain as composition sugar or as a specialized pipeline form

#### Environment model in `lower.py`
Use an explicit lexical environment for lowering:
- map workflow identifiers to IR references or imported function nodes
- distinguish imported names from local bound names
- keep shadowing rules explicit

Suggested approach:
- imported names enter the environment first
- later `Bind`s extend the environment lexically
- `Name("x")` refers to the nearest lexical binding

#### Non-goals for the first version of `lower.py`
Do not implement yet:
- execution DAG creation
- interpolation resolution
- bash validation
- forcing/saturation analysis
- closure simplification
- aggressive desugaring of chain into lambdas

The goal of the first lowering pass is only to produce a faithful, lazy semantic IR.

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

## Next implementation steps

1. Tighten workflow semantic checking around function-valued workflows further:
   - reduce remaining approximation gaps between lambda inference and top-level workflow signature inference
   - make lazy partial-application modeling less field-set-oriented and more value/provenance-oriented
   - decide how closures should carry bound values/provenance, not just bound field names

2. Normalize semantic error reporting:
   - move toward a single `errors` surface instead of legacy `chain_errors` naming
   - keep diagnostics output aligned with that simplification

3. Continue semantic IR cleanup:
   - decide whether `Closure` remains its own node or is folded into a function-like representation later
   - consider whether repeated `Apply(Function, arg)` structure should be shared more aggressively before forcing

4. Implement `python/swl/ir/force.py`:
   - lower semantic IR into an execution DAG / graph form
   - instantiate cached workflow bodies only when forced
   - realize task applications as execution leaves
   - preserve lexical environment and lazy application provenance while forcing

### Planned design for `python/swl/ir/force.py`

#### Purpose
`force.py` should be the boundary between:
- semantic IR as a lazy tree of values/functions/applications
- execution IR as a realized graph of concrete computations and data dependencies

It should not redo parsing or semantic validation.
Its job is to:
- force a chosen semantic IR root
- instantiate lazy workflow/task applications only when needed
- build a graph/DAG with sharing and explicit dependencies

#### Inputs and outputs
Suggested input:
- a semantic IR node from `python/swl/ir/lower.py`
- optionally an initial environment of bound values
- optionally a cache/registry of already-forced subcomputations

Suggested output:
- a forced execution graph rooted at one or more output value nodes
- enough metadata to print/debug the graph
- stable node identities for repeated/shared subcomputations
- a self-contained expanded DAG template with imported tasks/workflows fully resolved
- explicit external input placeholders for user-supplied runtime data

#### Core forcing rule
Forcing should be demand-driven.
It should only realize nodes needed to obtain the demanded result.

High-level cases:
- `Literal`, `Record`, `Field`, `Update`: force structurally / recursively as values
- `Function`, `Lambda`: remain values until applied and forced
- `Apply`: do not execute immediately; inspect callee and saturation state first
- `Chain`: realize as composed applications/connections between callable values
- `Block`: evaluate bindings lazily in lexical environment, then force the final result

#### Environment model during forcing
Use an explicit lexical environment separate from IR nodes.
Recommended shape:
- `ForceEnv(parent=None, values={})`

Responsibilities:
- lexical name lookup for `Name`
- lazy binding of `Bind` values
- shadowing according to semantic scope rules
- allow workflow body instantiation with captured/bound values

#### Execution graph node families
A first-pass execution IR can stay small.
Suggested node families:

- `ExecLiteral(value)`
- `ExecRecord(fields)`
- `ExecField(source, name)`
- `ExecMerge(left, right)`
- `ExecTaskCall(function, arg, outputs)`
- `ExecValueRef(node, field=None)`

Possible later additions:
- `ExecWorkflowInstance(...)`
- `ExecClosure(...)`
- explicit dependency edge objects if a purely node-based structure becomes awkward

The important property is graph identity and sharing, not exact class names yet.

#### Forcing task applications
When forcing an `Apply(Function(kind='task'), arg)`:
1. normalize the argument to the task input model
2. perform scalar-to-record lifting if needed
3. merge with any already-bound partial inputs
4. if still unsaturated, keep a function/closure-like forced value rather than emitting execution
5. if saturated, emit one fully-resolved `ExecTaskCall` planning node
6. expose task outputs as graph values/refs rather than executing bash immediately

For the compiled DAG artifact:
- task imports should no longer remain as unresolved imports
- each task call node should carry the task identity and enough task metadata for later execution without consulting source files again

This is where task interpolation defaults and runtime parameter defaults will later become relevant.

#### Forcing workflow applications
When forcing an `Apply(Function(kind='workflow'), arg)`:
1. normalize and merge the argument with any bound inputs
2. if unsaturated, return a closure/function-like value
3. if saturated, instantiate the cached workflow body in a fresh lexical environment
4. bind the workflow parameter to the merged argument value
5. force the instantiated body result
6. memoize by workflow function identity + bound input identity when possible

Key point:
- forcing a workflow should inline/instantiate its body only at force time, not earlier
- the final compiled DAG should be self-contained, so imported workflow boundaries should not remain as unresolved references in the serialized result

#### Forcing lambda applications
When forcing an `Apply(Lambda(...), arg)`:
- create a fresh lexical environment
- bind the lambda parameter to the argument value
- force the lambda body in that environment

This should follow the same mechanics as forcing an imported workflow body, except the body is already present locally.

#### Partial application representation during forcing
Forcing still needs a runtime-level representation for unsaturated callables.
A first implementation can use a small helper value like:
- `ForcedFunction(function, bound_arg=None)`

Responsibilities:
- remember original callable
- remember already-bound argument structure/value
- support another round of application before any execution node is emitted

This keeps forcing faithful to lazy partial application.

#### Sharing and memoization
`force.py` should introduce graph sharing aggressively where safe.
At minimum, memoize:
- forced imported workflow bodies by function identity
- saturated task calls by function identity + normalized bound argument identity
- projections from the same produced record when practical

This avoids DAG blow-up and is one of the main reasons `force.py` should exist separately from semantic lowering.
Sharing/memoization is an internal compilation optimization only; the final serialized DAG should still be self-contained and not depend on source-file reloading.

#### Chain forcing
`Chain(items)` should force as composition, not as immediate execution of every item in isolation.
A practical first rule:
- treat a chain as left-to-right composition of callables
- forcing a chain as a value yields a callable/composed callable
- forcing a saturated application of that chain realizes the underlying task/workflow calls in order

If easier initially, `force.py` may first desugar `Chain` into nested workflow-style application/composition before graph construction.

#### Field projection and record updates
`Field` and `Update` should remain lazy/value-level until they depend on concrete forced computations.
Examples:
- projecting `bam` from a saturated task call should produce a reference to that output of the `ExecTaskCall`
- updating two forced records should produce a merged value view, not duplicate execution

#### Serialization / compiled artifact
The forced DAG should be suitable for writing to disk as a compiled artifact.
Requirements:
- self-contained: all imported workflows expanded and all task references resolved into concrete task-planning nodes
- explicit external input placeholders for user-supplied runtime data
- no dependency on source-file lookup at execution time
- deterministic node/edge structure suitable for stable serialization

Recommended first format:
- JSON, for readability and testability

Possible later formats:
- msgpack
- protobuf
- custom binary / bytecode-like encoding

#### Proposed module structure
Initial contents of `python/swl/ir/force.py` should likely include:
- small forced-value helper classes
- small execution-node dataclasses
- `Forcer` class with caches and lexical environment helpers
- `force(node)` / `force_file(path)` entrypoints
- helpers to produce a serializable DAG form
- optional pretty-printer or separate `eval_force.py` later

Current implemented forcing shape:
- root callable workflows, including root chain-valued workflows, are auto-instantiated with symbolic external inputs
- compiled JSON is executor-oriented rather than syntax-oriented
- serialized payload keeps only minimally required executor data:
  - external `inputs` keyed by name
  - per-task `id`, `name`, `path`, `deps`
  - per-task normalized input bindings
  - per-task normalized declared outputs
  - per-task normalized run params
  - per-task embedded `script`
  - final workflow `outputs`
- raw task annotation sections are intentionally not serialized
- dependency traversal is available directly via each task's `deps`
- source-file rereads should not be required at execution time

Remaining force work from here:
- improve workflow-input metadata so compiled external inputs retain more precise type/description information in more cases
- finalize executor contract for evaluating interpolation/default expressions against runtime bindings
- decide whether compiled task `path` remains required provenance only or also part of executor identity
- canonicalize final output formation so workflows written with explicit task applications and record updates compile to the same output shape as equivalent chain sugar

Specific canonicalization target discovered during compilation:
- `tests/pipe.swl` and `tests/function.swl` are intended to be equivalent
- today `pipe.swl` lowers to `ir.Chain(...)`, while `function.swl` lowers to explicit `Apply(...)` + `Update(...)` structure
- the current forcing pass already knows how to merge chain outputs into a flat output map, but it does not yet flatten general nested record merges at the workflow boundary
- forcing should therefore flatten final `Merge(Record(...), ...)` structures into a canonical flat output map using right-biased update semantics
- after that, both chain sugar and explicit update-based workflow composition should serialize the same final outputs

Possible later lowering cleanup in `python/swl/ir/lower.py`:
- add a small canonicalization pass after lowering expressions but before forcing
- detect explicit workflow-composition patterns such as:
  - `a = f x`
  - `b = g (x // a)`
  - final result built from `a // b` (and longer left-to-right repetitions)
- rewrite those patterns into a canonical chain-like representation, or into a normalized record-producing block shape
- this is optional once `force.py` flattens final merges correctly, but it would reduce divergence between equivalent source forms earlier in the pipeline and make equality/testing/debugging simpler

## Current responsibility split: `lower.py` vs `force.py`

### `python/swl/ir/lower.py`
`lower.py` should be the semantic normalization boundary.
Its job is to erase parser-specific syntax details and produce a smaller, more canonical semantic IR.
That means it should own:
- import normalization (`import` syntax -> callable/function IR values)
- lexical block/lambda/binding normalization
- canonicalization of equivalent source forms where no execution forcing is required
- attachment/preservation of signature metadata where available
- structural simplifications that are semantics-preserving and do not depend on runtime data or DAG identity

In other words, `lower.py` should answer:
- what is the program structurally?
- which forms are semantically equivalent before any execution planning happens?

### `python/swl/ir/force.py`
`force.py` should be the execution-planning boundary.
Its job is to take the semantic IR and instantiate it into an executor-facing DAG template.
That means it should own:
- symbolic external input creation for root callable workflows
- saturation / partial-application handling at compile time
- imported workflow expansion at force time
- concrete task-call node creation
- task-call memoization / graph sharing
- dependency extraction and final executor-facing JSON serialization

In other words, `force.py` should answer:
- given the normalized semantic meaning, what concrete execution template should an executor traverse?

### Problem with the current split
Right now `force.py` is compensating for too much structural divergence that should ideally have been reduced earlier.
`pipe.swl` and `function.swl` are the clearest example:
- `pipe.swl` lowers directly to a chain-oriented semantic form
- `function.swl` lowers to a block of explicit applications plus nested record updates
- both are intended to mean the same workflow composition
- because the lowerer preserves those different shapes, forcing has to reconstruct their equivalence late

This is why forcing currently carries logic for:
- merging chain outputs
- flattening final output records/merges
- recovering input metadata from emitted tasks
- composing signatures that the lowerer could potentially make more explicit earlier

## How to get `lower.py` to do more work

The right direction is not to make `lower.py` eagerly build DAG nodes.
Instead, it should do more semantic canonicalization so `force.py` sees fewer distinct-but-equivalent shapes.

### 1. Add an explicit normalization pass inside `lower.py`
Recommended structure:
- first lower parser AST to semantic IR
- then run `normalize(node)` before returning

This keeps the parser-to-IR translation simple while giving a dedicated place for semantic canonicalization.
That pass can work recursively over the semantic IR and rewrite equivalent trees to a more stable normal form.

### 2. Normalize record-update structure more aggressively
Beyond closed-record merging, the lowerer should normalize update trees into a predictable shape.
Examples:
- associate nested updates consistently
- flatten statically-known record fragments where safe
- normalize `Update(Update(a, b), c)` into one ordered update chain representation, even when some pieces are not closed records

The goal is not to erase open-record semantics, but to remove purely syntactic nesting noise.
This would make final-result assembly in explicit workflows much more uniform before forcing.

### 3. Recognize explicit stage-building patterns in blocks
Without converting everything into `ir.Chain`, the lowerer can still recognize a canonical staged-composition shape.
For workflows like:
- `a = f x`
- `s = g (x // a)`
- `c = h (x // a // s)`
- final result `a // s // c`

the lowerer can attach or produce a normalized staged form that records:
- stage order (`a`, `s`, `c`)
- each stage's callable
- each stage's argument environment as a progressive accumulation of previous stage outputs plus root input
- final output assembly as the union of stage result records

Important point:
- this does **not** have to rewrite the program into `ir.Chain`
- instead it can canonicalize it into a normalized block/result shape that is equivalent to chain sugar but still preserves the explicit style

A practical way to do this is:
- detect a block whose bindings are stage applications
- verify that each stage argument is the root input updated with previously-bound stage results
- verify that the final result is an update-union of those stage bindings
- then rewrite just the *result assembly* and stage metadata into a canonical staged normal form

### 4. Preserve progressive workflow-composition metadata in semantic IR
If full structural rewriting feels too invasive, the lowerer can still enrich the semantic IR with canonical metadata.
Examples:
- attach a normalized composed signature to staged blocks
- annotate canonical stage order for certain normalized blocks
- preserve explicit information about which updates are stage-result unions rather than arbitrary record merges

This would allow `force.py` to consume a more informative semantic IR without having to rediscover everything structurally.

### 5. Normalize final result construction
One concrete target is the workflow final expression.
Today:
- chain sugar tends to arrive at force time already looking like a union of stage outputs
- explicit composition arrives as nested updates/merges

The lowerer should normalize final workflow result expressions so that equivalent programs expose the same logical output assembly before forcing.
That means:
- a final union of stage result records should lower to one canonical form regardless of whether it came from syntactic chain sugar or explicit `//` expressions
- forcing should then only have to serialize that canonical result, not rediscover its meaning

## Desired end state for `tests/function.swl` vs `tests/pipe.swl`

The goal is:
- same task set
- same task dependencies
- same external inputs
- same workflow outputs
- almost identical compiled JSON

The remaining acceptable differences should be limited to provenance/debug details, not semantic structure.
Specifically, we do **not** want to merely special-case `function.swl` by converting it into `ir.Chain`.
Instead, we want a more general canonicalization story where:
- chain sugar
- explicit stage bindings
- explicit final result unions

all lower into one semantically stable normal form that `force.py` can plan from uniformly.

## Proposed implementation plan

### Phase A: strengthen normalization in `lower.py`
1. Add `normalize(node)` after lowering.
2. Normalize update trees into a stable internal shape.
3. Detect stage-building workflow blocks/lambdas that represent progressive composition.
4. Canonicalize their final result assembly into a stable union-of-stage-results form.
5. Preserve or attach normalized composition metadata/signatures where recoverable.

## Redesign direction: identity-based semantic IR for DAG-friendly lowering

The current `Compose`-based normalization helped align `pipe.swl` and `function.swl`, but it still carries a tree-oriented weakness: later stage arguments rebuild nested `Update(...)` structure instead of referring explicitly to previously defined semantic values.
That is acceptable as a transitional normalization, but it is not the best long-term design for scope-correct DAG creation.

If we are not constrained to keep changes minimal, the better design is to make semantic identity explicit in the IR.
The real problem is not just pretty normalization of update trees; it is that semantic sharing is currently implicit and tree-shaped, so `force.py` must recover graph structure later.
The redesign should move that sharing earlier into the IR itself.

### Core principle
Represent bound computations once, with stable identity, and refer to them explicitly everywhere else.
Do **not** substitute full bound value trees into later expressions.
Do **not** rely on structural equality of trees to recover sameness.
Instead:
- each semantic binding gets a unique identity
- later uses become explicit references to that identity
- `force.py` memoizes by semantic identity rather than by syntactic shape

## Proposed semantic IR redesign

### New binding/reference nodes
Add explicit identity-bearing nodes to `python/swl/ir/node.py`.

#### `Variable`
Represents one bound semantic computation.
Suggested fields:
- `id: int`
- `name: str`
- `value: Node`

This is the semantic definition site of a computed value.
Examples:
- `a = align x`
- `s = sort (x // a)`
- `c = call (x // a // s)`

should each become one `Variable`.

A workflow binding should lower directly to a `Variable`; there should not be a separate semantic `Bind` layer once lowering is complete.

#### `Ref`
Represents a use of a previously defined semantic value.
Suggested fields:
- `id: int`
- `name: str`

This is the semantic use site of a computed value.
A later stage argument should refer to prior stage results through `Ref`, not by copying the original `Apply(...)` tree.

### Updated block structure
`Block` should hold identity-bearing variables rather than anonymous binds whose values get substituted later.
That means `Block` should conceptually become:
- a sequence of `Variable`s
- plus a result expression that may reference them through `Ref`

The parser AST may still have bind syntax, but semantic lowering should produce `Variable` directly.

## Lowering strategy under the redesign

### Current problem in `lower.py`
Today lowering substitutes environment values directly into later expressions.
That means later expressions contain duplicated copies of earlier computations.
This is the root cause of tree duplication and ambiguous sharing.

### New lowering rule
When lowering a binding:
1. allocate a fresh semantic variable id
2. lower the bound expression once
3. emit `Variable(id, name, value)`
4. place `Ref(id, name)` in the environment for later uses

Then:
- the binding is defined exactly once
- later references preserve scope-correct identity
- nested scopes are safe because they allocate fresh ids even for shadowed names

### Example target shape
For:

```swl
\x ->
    a = align x
    s = sort (x // a)
    c = call (x // a // s)
    a // s // c
```

the lowered semantic IR should conceptually look like:

```python
Lambda(
  param='x',
  body=Block(
    bindings=[
      Variable(id=1, name='a', value=Apply(Function('align'), Name('x'))),
      Variable(id=2, name='s', value=Apply(Function('sort'), Update(Name('x'), Ref(1, 'a')))),
      Variable(id=3, name='c', value=Apply(Function('call'), Update(Update(Name('x'), Ref(1, 'a')), Ref(2, 's')))),
    ],
    result=Update(Update(Ref(1, 'a'), Ref(2, 's')), Ref(3, 'c')),
  ),
)
```

The exact syntax can vary, but the essential property is:
- definitions occur once
- uses are explicit references
- no repeated `Apply(align, x)` subtrees are copied into later arguments

## Relationship to `Compose`

### Transitional view
`Compose` is still useful as a normalization concept because it captures staged workflow composition explicitly.
However, in the long-term design it should be built on top of value identity, not on top of copied update trees.

### Better `Compose` model
If `Compose` remains, its stages should refer to semantic values explicitly.
That means a stage should represent a stage result variable, not merely a callee/arg pair with ambiguous meaning.
A cleaner stage model would be something like:
- a variable/result identity
- the callee
- the argument expression, which may include `Ref`s

or even simply a `Variable` specialized to composition staging.

### Important conclusion
The real long-term goal is **not** merely to make `Compose` prettier.
The goal is to make the entire lowered semantic IR explicitly shareable and scope-correct.
Once that exists, `Compose` becomes optional or at least much simpler.

## Forcing strategy under the redesign

### Current limitation in `force.py`
`force.py` currently reconstructs DAG sharing from tree structure, task-input normalization, and ad hoc memoization.
That is workable but indirect.

### New forcing rule
`force.py` should force semantic definitions lazily by identity.
That means:
- maintain a mapping from `Variable.id` to forced result
- when a `Ref(id, ...)` is encountered:
  - if already forced, reuse the existing forced value
  - otherwise force the referenced definition once and cache it

This turns forcing into graph evaluation over semantic let-bindings rather than repeated tree traversal.

### Benefits
- no duplicated task emission for repeated semantic stage uses
- scope correctness by construction
- much easier DAG memoization
- cleaner equivalence between syntactic sugar and explicit workflow definitions

## Why this is better than structural caching

### Structural caching is insufficient
Caching by raw lowered tree shape like:
- `Update(Name('x'), Name('a'))`
- or by repeated `Update(...)` object interning

is not enough because:
- it is not inherently scope-safe
- it conflates structural similarity with semantic identity
- it still relies on tree reconstruction heuristics

### Identity-based refs are scope-safe
With `Variable` / `Ref`:
- two shadowed names in nested scopes get different ids
- repeated uses of the same binding share one semantic identity
- memoization can safely key on that identity

This is the correct abstraction boundary for DAG creation.

## Implications for `pipe.swl` vs `function.swl`
The ideal end state is that both forms normalize to the same identity-based semantic graph.

### `pipe.swl`
Chain syntax should normalize into explicit staged semantic definitions/refs.
Conceptually, `align | sort | call` should normalize to the same staged let-graph as the explicit form.

### `function.swl`
The explicit function form already expresses those stages directly; lowering should preserve them as explicit semantic defs/refs rather than flattening them into copied trees.

Then both workflows converge before forcing, and the forced DAG becomes naturally identical.
This is better than special-casing chains or relying on forcing-time reconstruction.

## Recommended redesign phases

### Phase 1: introduce semantic identity
1. Add `Variable` and `Ref` to `python/swl/ir/node.py`.
2. Update `lower.py` so binding lowering allocates unique ids.
3. Lower bindings directly to `Variable`.
4. Stop substituting bound value trees into later expressions; store `Ref`s in the env instead.
5. Update repr/debug printing accordingly.

### Phase 2: normalize composition on top of refs
1. Revisit `Compose` so stages refer to identity-bearing values.
2. Normalize both `Chain(...)` and explicit staged blocks into the same ref-based staged representation.
3. Keep update normalization only as a helper, not as the primary sharing mechanism.

### Phase 3: rewrite forcing around semantic identity
1. Teach `force.py` to resolve `Ref` through variable definition lookup.
2. Memoize forced results by semantic id.
3. Emit each semantic application/task at most once.
4. Reduce forcing heuristics that currently reconstruct sharing from tree shape.

### Phase 4: simplify after convergence
Once the ref-based design is in place:
- simplify `Compose` if it becomes redundant
- simplify chain-specific forcing logic
- simplify update-flattening logic used only for equivalence reconstruction

## Target end state
The semantic IR should be a tree-with-refs or graph-shaped representation of workflow meaning, not a pure copied tree.
`force.py` should then be a mostly straightforward graph-to-DAG pass.

That is the best long-term route to:
- no duplicated DAG work
- scope-correct caching/memoization
- clean equivalence between syntactic sugar and explicit workflows
- less forcing-time guesswork

### Phase B: simplify `force.py` once lowering is stronger
After stronger lowering normalization, `force.py` should be able to rely more on the lowered shape and do less equivalence reconstruction.
That should allow simplifications in:
- final output flattening logic
- chain-vs-explicit-composition special handling
- signature composition heuristics

### Phase C: regression tests
Add tests asserting that compiling:
- `tests/pipe.swl`
- `tests/function.swl`

produces nearly identical JSON, especially for:
- `inputs`
- `tasks`
- `outputs`

If exact object equality is too strict, compare normalized serialized forms with provenance-only differences ignored.

#### Suggested implementation order
1. Define tiny execution-node dataclasses.
2. Force literals/records/fields/updates/names/blocks.
3. Add callable forcing for `Function` and `Lambda` as non-executing values.
4. Add partial-application handling.
5. Add saturated task-call forcing to graph nodes.
6. Add saturated workflow/lambda instantiation.
7. Add chain forcing/composition.
8. Add memoization and graph-printing.

#### Initial tests for `force.py`
Add focused unit tests covering:
- forcing a saturated imported task application produces one task-call node
- forcing a partial task application does not produce a task-call node yet
- forcing a saturated imported workflow application instantiates its cached body
- forcing repeated use of the same imported function shares cached graph/function objects where intended
- forcing field access on task outputs produces output references, not duplicated task calls
- forcing a chain realizes underlying calls in left-to-right dependency order
- serialized compiled DAG contains no unresolved workflow/task import references
- serialized compiled DAG can be loaded without reading original source files again

5. Add more semantic tests around scope and function-ness:
   - nested block duplicate-binding errors
   - uniform shadowing behavior across nested scopes
   - imported workflow final-value failures
