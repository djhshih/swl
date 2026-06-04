# SWL Architecture

SWL is a workflow language that compiles to a logical DAG of task invocations. The pipeline has these stages:

1. **Workflow Lexing** → token stream
2. **Workflow Parsing** → AST (`syntax/wf/node.py`)
3. **Task Parsing** → parsed task annotation + bash body (`syntax/task/node.py`)
4. **Semantic Checking** → validated imports, scope, types, inferred signature (`semantic/wf/check.py`)
5. **IR Lowering** → semantic IR tree (`ir/node.py`)
6. **IR Forcing** → concrete DAG (`dag/`)
7. **Output** → JSON-serialized DAG emission

---

## Package Structure

The compiler lives under `python/swl/` as a single package:

```
swl/
├── api.py              # Public API (compile, force, load, transpile)
├── compile.py          # CLI entry point
├── loader.py           # File cache, parsed-task & workflow memoization
├── repl.py             # Interactive REPL
├── types.py            # Type normalization & CWL/WDL/NF type mapping
├── syntax/wf/          # Workflow language lexing/parsing
├── syntax/task/        # Task script annotation/body parsing
├── semantic/wf/        # Workflow semantic analysis
├── semantic/task/      # Task type system
├── ir/                 # Intermediate representation (nodes + lowering)
├── dag/                # DAG construction, forcing, finalization
├── eval/               # Debug CLI tools (6 scripts)
└── transpile/          # Backend code generators (CWL, WDL, Nextflow)
```

---

## Stage 1: Workflow Lexing

**File:** `syntax/wf/lexer.py`

**Input:** Raw SWL source string (`.swl` file).

**Output:** A stream of `Token` objects with a `TokenType` enum.

**What happens:**
- `\r\n` → `\n`, `\r` → removed (per spec).
- Indentation is tracked via an `indent_stack`. Each block of 4 spaces (or a tab) at the start of a line emits a `bstart` token when indentation increases, and `bend` tokens (one per level) when it decreases. This converts the indentation-based block syntax into explicit bracketing tokens.
- Comments (`#` to end of line) are discarded.
- Non-newline whitespace is consumed silently (the spec says "removed" but the lexer achieves this by skipping it).
- Multi-character tokens (`->`, `//`) are recognized before single-character fallbacks.
- String literals (`"..."`) handle escaped quotes.
- When `{` or `|` or `,` or `->` is encountered, `ignore_eol` is set so the next newline is suppressed (these constructs are inline).

**Token types produced:** `id`, `str`, `num`, `equal`, `colon`, `comma`, `dot`, `update` (`//`), `chain` (`|`), `lparen`, `rparen`, `lbrace`, `rbrace`, `bslash` (`\`), `arrow` (`->`), `bstart`, `bend`, `eol`, `eof`.

**Representation at this stage:** A flat list of `Token` objects, e.g. `[bslash, id("x"), arrow, bstart, id("y"), equal, id("f"), id("x"), bend, eof]`.

---

## Stage 2: Workflow Parsing

**Files:** `syntax/wf/parser.py`, `syntax/wf/node.py`

**Input:** Source string (re-lexed internally via a `Lexer` instance).

**Output:** An AST of `Expr` subclass instances. Each `Expr` has a `type` field from `NodeType` enum.

**Node types in the AST:**
| NodeType | Class    | Meaning |
|----------|----------|---------|
| `id`     | `Identifier` | Variable reference |
| `str`    | `String`     | String literal |
| `num`    | `Number`     | Number literal |
| `block`  | `Block`      | Sequence of expressions, each separated by newline |
| `bind`   | `Binding`    | `name = expr` |
| `rec`    | `Record`     | `{ key: expr, ... }` |
| `get`    | `Get`        | `expr.member` (record field access) |
| `update` | `Update`     | `left // right` (record merge) |
| `fun`    | `Function`   | `\param -> body` (lambda) |
| `apply`  | `Apply`      | `fun arg` (function application) |
| `chain`  | `Chain`      | `left | right` (pipeline) |

**What happens:**
- The parser is a recursive descent parser that follows the GBNF grammar from `spec.md`.
- Precedence is handled by nesting parse methods: `_parse_chain_expr` (lowest) → `_parse_update_expr` → `_parse_apply_expr` → `_parse_get_expr` → `_parse_term` (highest).
- A block is parsed via `_parse_block()`, which collects expressions separated by `eol` tokens. It validates that the final expression is not a binding.
- Indented sub-blocks (lambda bodies) are parsed by `_parse_block(inner=True)`, which expects `bstart` and `bend`.
- `_parse_function()` handles `\` name `->` body.
- `_parse_record()` handles `{ key: value, ... }` with comma separators.
- The `import` call is detected later (not at parse time) — at parse time it's just an `Apply` of `id("import")` to a `String`.
- **Shared pattern-matching** (`syntax/wf/builtins.py`): Functions `match_import`, `match_map`, and `match_map_by` are shared between the checker and lowerer to detect `import("path")`, `map f xs`, and `map_by f key xs` AST patterns.

**Representation at this stage:** A tree of `Expr` objects rooted at a `Block`. Example for `\x ->\n    y = f x\n    y`:
```
Block([
  Function(
    param=Identifier("x"),
    body=Block([
      Binding(Identifier("y"), Apply(Identifier("f"), Identifier("x"))),
      Identifier("y")
    ])
  )
])
```

---

## Stage 3: Task Parsing

**Files:** `syntax/task/parser.py`, `syntax/task/node.py`

**Input:** A `.sh` file containing comment-based annotations followed by a bash script body.

**Output:** A `Task` object with `.annotation` (an `Annotation` with `.doc` string and `.sections` list) and `.body` (the raw bash script text).

**What happens:**
- Lines are split; leading `#` lines are annotation, the rest is the script body.
- `#` prefix is stripped along with one optional space.
- Annotation structure:
  - First non-blank annotation line must start with `@` — this is the doc string.
  - Then `in`, `out`, or `run` sections, each containing parameter lines.
  - Each parameter line has: comma-separated names, optional type (`file`, `str`, `int`, `float`, optionally with `?` or `[...]`), optional default (`= value`), optional description (`| text`).
  - Default values are parsed by the interpolation parser.
  - Description continuations start with `|` on the next line.

**Supporting file — `syntax/task/interpolation.py`:**
Parses shell-style string interpolation (`${var}`, `$var`, `${expr}`) into a `Word` containing `Literal`, `Var`, or `Expr` parts. Used for default values in task parameters and for variable mentions in bash commands.

**Supporting file — `syntax/task/bash.py`:**
Parses the bash body into a `Script` of `Assignment` and `Command` statements. Each `Command` records which words contain interpolations. Used by the CWL transpiler to identify variables referenced in the script.

**Representation at this stage:** A `Task` object:
```
Task(
  annotation=Annotation(
    doc="align reads",
    sections=[
      Section(kind=IN, params=[Param(names=["fastq1","fastq2"], type="file", default=..., desc=...)]),
      Section(kind=OUT, params=[Param(names=["bam"], type="file", default=Word(...), desc=...)]),
      Section(kind=RUN, params=[Param(names=["cpu"], type="int", default=Word([Literal("2")]))]),
    ]
  ),
  body="bwa mem ...\n"
)
```

---

## Stage 4: Semantic Checking — Task Types

**File:** `semantic/task/type.py`

**Input:** A parsed `Task` (from stage 3).

**Output:** A `TaskSignature` with typed `Param` dictionaries for `inputs`, `outputs`, and `run`.

**What happens:**
- `signature_from_task()` walks the annotation sections and builds `Param` objects with resolved `TypeKind` enums.
- `TypeKind` values: `FILE`, `STR`, `INT`, `FLOAT`, `MEMORY`, `TIME`, and optional/array variants.
- Run parameters (`memory`, `time`, `cpu`, `image`) get special normalization:
  - `memory` is parsed as a literal like `8G` → integer MB.
  - `time` as `HH:MM:SS` or `D-HH:MM:SS` → integer minutes.
  - `cpu` as integer.
  - `image` as string.
- `TypeChecker` class provides `check_chain(left, right)` which validates that common-named outputs of `left` are type-compatible with inputs of `right`. The compatibility matrix follows the spec (e.g., `file` → `file?` is OK, but `file` → `str` is not).

---

## Stage 5: Semantic Checking — Workflow

**Files:** `semantic/wf/check.py`, `semantic/wf/imports.py`, `semantic/wf/scope.py`, `semantic/wf/infer.py`, `semantic/wf/signature.py`, `semantic/wf/bashvars.py`

**Input:** Parsed workflow AST (from stage 2) plus imported tasks/workflows.

**Output:** A `WorkflowCheck` object containing: validated AST, imports, errors, inferred inputs, workflow signature, workflow type, batch-vs-simple classification.

**What happens:**

1. **Import resolution** (`semantic/wf/imports.py`): Scans top-level bindings for `name = import("path")` calls. `.sh` paths → task imports; `.swl` paths → workflow imports (recursively loaded). Circular imports are detected via a `_loading` stack on the `Loader` (`loader.py`).

2. **Scope checking** (`semantic/wf/scope.py`): Walks the AST ensuring no duplicate variable names within the same scope. Also checks for forward references (a binding referencing a later binding in the same block) and self-references (a binding referencing itself). Nested scopes (lambda params, inner blocks) may shadow outer names. Chain type checking (`check_chains`) calls the task-level `TypeChecker.check_chain()` for known tasks in `A | B` pipelines.

3. **Bash variable validation** (`semantic/wf/bashvars.py`): Validates that all variable references in imported task bash bodies are resolvable from declared inputs or run parameters.

4. **Input inference** (`semantic/wf/infer.py`): Evaluates the workflow body in a symbolic/semantic fashion to determine what inputs the workflow demands. Uses a set of semantic value types:
   - `OpenRecord` / `ClosedRecord` — record values with known/unknown fields.
   - `FunctionValue` — a known function with a signature.
   - `ClosureValue` — a partially applied function.
   - `ComputationValue` — a fully applied function (a task invocation).
   - `TableValue` — a table (columnar batch data).
   - `UnknownValue` — unresolved.
   - `TypedValue` — a value with a known type but unknown origin.

   The semantic evaluator (`_eval_expr`) recursively walks the AST. When it encounters a field access like `x.field` on an `OpenRecord`, it adds `field` to the `demanded` set (this is how input fields are inferred). When it encounters a function application, it checks saturation — if all required inputs are available, it produces a `ComputationValue`; otherwise a `ClosureValue`.

5. **Signature construction** (`semantic/wf/signature.py`): Determines whether the workflow is batch (`tab -> tab`) or simple (`rec -> rec`) based on whether `map`/`map_by` is used. Builds the final `TaskSignature` from inferred inputs and outputs. Also produces `workflow_type` as a `FunctionType(input_type, output_type)` from `semantic/wf/type.py`.

**Type AST** (`semantic/wf/type.py`):
- `ScalarType("file" | "str" | "int" | "float" | "?")`
- `ArrayType(item)`
- `RecordType(fields, open)`
- `TableType(columns)`
- `FunctionType(input, output)`

**Loader** (`loader.py`): The `Loader` class provides in-memory file caching, parsed-task memoization (`_parsed_tasks`), and checked-workflow memoization (`_checked_workflows`). All file reading goes through `read_file()` which supports an in-memory `files` dict (used by tests) as well as real filesystem access. The `Checker` and `Lowerer` both hold a `Loader` instance.

---

## Stage 6: IR Lowering

**Files:** `ir/lower.py`, `ir/node.py`

**Input:** Workflow AST + imports + signature from the semantic checker.

**Output:** A semantic IR tree of `ir.Node` dataclass instances (all frozen/immutable).

**IR node types:**
| Node       | Meaning |
|------------|---------|
| `Literal`  | A literal value (string, number) |
| `Input`    | A workflow input placeholder |
| `Name`     | An unresolved name reference (will become a workflow input or local variable) |
| `Ref`      | A resolved reference to a `Variable` by id |
| `Record`   | A record literal `{ field: value, ... }` |
| `Field`    | Field access `record.name` |
| `Update`   | Record merge `left // right` |
| `Function` | An imported task/workflow with name, kind, signature, path, optional body |
| `Lambda`   | A lambda with param name, body Block, optional signature |
| `Closure`  | A function with a partially bound argument |
| `Apply`    | Function application |
| `Map`      | Special map/map_by node |
| `Variable` | A let-binding: `id = value` |
| `Block`    | A sequence of bindings and a final result expression |
| `Unknown`  | Placeholder (should be eliminated during lowering) |

**What happens:**
- The `Lowerer` uses `Checker.load()` to semantically check the workflow, then walks the AST to produce IR.
- **Imports** become `ir.Function` nodes carrying their kind (`task`/`workflow`), path, and signature. Workflow imports recursively lower their body via `_cached_workflow_body()`.
- **Lambdas** become `ir.Lambda` with a `Block` body.
- **Bindings** become `ir.Variable` nodes with unique integer `id`s. Subsequent references are `ir.Ref` nodes pointing to the `Variable` by id.
- **`map`/`map_by` detection**: The lowerer uses `builtins.match_map` / `builtins.match_map_by` (from `syntax/wf/builtins.py`) to detect `map f xs` and `map_by f key xs` (which are `Apply(Apply(id("map"), f), xs)` in AST) and collapses them into `ir.Map` nodes with an optional `key` field.
- **Pipeline chains** (`A | B | C`) are desugared into a lambda block following the spec:
  ```
  \_input ->
      _s1 = A _input
      _s2 = B (_input // _s1)
      _s3 = C (_input // _s1 // _s2)
      _s1 // _s2 // _s3
  ```
- **Normalization** (`normalize()`): Recursively transforms the IR:
  - Lambda bodies are ensured to be `Block` nodes.
  - `Map` callables are "materialized" — lambdas and partial applies used with `map`/`map_by` are promoted to synthetic `Function` nodes with generated names and inferred signatures.
  - Update chains (`a // b // c`) are flattened and reconstructed as nested `Update` nodes.
  - Variable references are resolved through their binding maps.

**Representation at this stage:** A tree of `ir.Node` instances. For example, `align | sort` desugared:
```
Lambda(
  param="_input",
  body=Block(
    bindings=[
      Variable(id=1, name="_s1", value=Apply(Function("align"), Name("_input"))),
      Variable(id=2, name="_s2", value=Apply(Function("sort"), Update(Name("_input"), Ref(1, "_s1")))),
    ],
    result=Update(Ref(1, "_s1"), Ref(2, "_s2"))
  )
)
```

---

## Stage 7: IR Forcing (DAG Construction)

**Package:** `dag/` — 8 files that collectively convert the semantic IR into a concrete, serializable DAG.

| File | Purpose |
|------|---------|
| `node.py` | DAG data types: `Input`, `Literal`, `Record`, `TableSource`, `Field`, `Merge`, `StepCall`, `Output`, `OutputSpec`, `ForcedFunction`, `DAG`. JSON serialization via `DAG.to_dict()` / `DAG.from_dict()`. |
| `context.py` | `ForceEnv` — simple scoped variable environment (parent pointer + values dict). |
| `binding.py` | Binding serialization: `binding_to_dict()` / `binding_from_dict()` for DAG values. |
| `evaluator.py` | Core forcing engine: `force_value()` recursively evaluates IR nodes into DAG values. |
| `forcer.py` | Top-level orchestration: `ForceState` dataclass, `force()` and `force_file()` entry points. |
| `emit.py` | Step call emission: `_emit_task_call()`, `_emit_workflow_call()`, `_emit_mapped_step()`. |
| `merge.py` | Merge canonicalization: collapsing, flattening, value-key computation for deduplication. |
| `finalize.py` | DAG finalization: input metadata refinement, step binding flattening, output spec building, unused input pruning, validation. |
| `tooldefs.py` | Tool definition construction: `_tool_definition()` returns full task definition dict, `_workflow_definition()` materializes sub-workflow DAGs, `_force_root()` applies final function to synthetic inputs. |

### What the Forcer does

1. **`force_value(node, env)`** — Recursively evaluates IR nodes in a `ForceEnv`:
   - `Literal` → `dag.Literal`
   - `Name` → either a `ForcedFunction` (if it's a known variable bound to a function) or an `Input` (workflow input)
   - `Record` → `dag.Record` with forced fields
   - `Field` → `dag.Field` with attempted projection (resolves record fields when possible, otherwise creates a `Field` chain)
   - `Update` → checks for unsupported operands (StepCall, Input+Record), rejects with clear error; otherwise produces `dag.Merge` or record merge
   - `Function` / `Lambda` → `ForcedFunction` wrapper
   - `Apply` → tries to apply the function to the argument
   - `Map` → maps function over table source
   - `Block` → evaluates in a local `ForceEnv`
   - `Variable` → forces the value and caches it

2. **Function application** (`_apply` in `evaluator.py`):
   - If the function is a built-in (`map`, `map_by`), handles partial application (binding `f`, then `key`/`xs`).
   - If the argument is a mapped/table source, delegates to `_apply_mapped`.
   - For tasks: checks saturation — if all required inputs are available in the bound record, emits a `StepCall` via `_emit_task_call`. Otherwise returns a `ForcedFunction` (partial application).
   - For workflows: similar saturation check, emits `StepCall` with `type='workflow'`.
   - For lambdas: if the bound argument is a `Record`, evaluates the lambda body directly.

3. **Step call emission** (`emit.py`):
   - Normalizes task inputs by projecting named fields from the bound argument.
   - Creates a `StepCall` with unique ID, dependencies, and tool definition.
   - Caches step calls by function path + input keys to deduplicate.

4. **Map emission** (`_emit_mapped_step` in `emit.py`):
   - Creates a `StepCall` with `map` metadata: source type (table/input), columns, grouping key for `map_by`.
   - Classifies step inputs into `scatter` (vary per table row) and `broadcast` (same for all rows). Inputs matching table columns are scattered; non-matching inputs are broadcasted. For generated workflows with a single record input, all inputs are scattered.
   - Validation: `DAG._validate_mapped_ports` ensures every input in `input_schema` belongs to exactly one of `scatter` / `broadcast`.

5. **Root forcing** (`_force_root` in `tooldefs.py`):
   - If the final value is a `ForcedFunction`, applies it to synthetic input placeholders to saturate it.
   - For batch workflows, creates a table input placeholder and applies `map`.

6. **DAG finalization** (`finalize.py`):
   - `_refine_input_metadata()`: Materializes workflow input metadata by merging inferred interface information with task input specs, including type, description, and optionality (detected via `?` suffix on type strings).
   - `_flatten_step_bindings()`: Expands `Record` and `Merge` bindings in step inputs into individual named bindings, flattening deeply nested merge trees via `_flatten_merge_value`.
   - `_final_outputs()`: Converts the forced root value into explicit top-level workflow outputs. `Merge` trees are canonicalized (record-record merges collapse); remaining merge trees are flattened.
   - `_build_output_specs()`: Wraps each output in an `OutputSpec` with inferred type, optionality, and the binding value.
   - `_assert_wireable_output()`: Rejects outputs that can't be expressed as `Input`, `Literal`, or `Field(Input|StepCall)`.
   - `_prune_unused_inputs()`: Removes input ports that aren't referenced by any step or output.
   - `dag.validate()`: Checks dependency DAG for cycles, unknown deps, and validates mapped-port scatter/broadcast classification.
   - Returns a `DAG(inputs, steps, outputs)` that satisfies the normalized DAG contract.

### Normalization obligations before final DAG emission

Before emitting the final DAG, the compiler must ensure all final-DAG invariants required by `dag.md` hold.

- **Semantic checking** is responsible for validating source-language meaning, inferring workflow inputs and outputs, and preserving source type information such as optionality in the inferred interface.
- **IR lowering** is responsible for desugaring source constructs into explicit IR forms, including pipeline expansion and explicit record/update structure.
- **IR forcing** is responsible for resolving applications into concrete step calls, saturating direct call interfaces, and carrying inferred interface information forward into concrete DAG values.
- **DAG finalization** is responsible for producing a fully materialized workflow interface and enforcing final normalized DAG invariants.

In particular, before DAG emission the compiler must:
- materialize workflow input and output metadata into explicit DAG interface specifications;
- preserve optionality from parsed task annotations and inferred workflow interfaces into DAG metadata;
- eliminate merge bindings from the final DAG by flattening them into explicit named bindings or outputs;
- saturate record literals that directly feed known task or workflow interfaces into named port bindings;
- ensure top-level workflow outputs are serialized in the normalized output form defined by `dag.md`; and
- reject any remaining binding form that falls outside the portable final DAG contract.

**Representation at this stage:** A `DAG` dataclass with concrete, serializable values. Each `StepCall` has:
- `id`: unique step name (e.g., `align`, `align_2`)
- `path`: path to the tool script
- `bindings`: dict of input name → normalized binding value
- `outputs`: list of output names
- `deps`: list of step IDs this step depends on
- `task`: the full tool definition (inputs, outputs, run, body, doc)
- `type`: `"task"` or `"workflow"`

The emitted DAG is the compiler artifact consumed by transpilers. Source-level constructs such as unresolved merges may exist in IR during lowering and forcing, but they must not survive into the final emitted DAG unless explicitly permitted by `dag.md`.

---

## Stage 8: Output

**File:** `compile.py` (CLI entry point), `api.py` (programmatic API), `dag/node.py` (serialization)

**Input:** A normalized `DAG` object.

**Output:** A JSON file written to `dag/<name>.json`.

**What happens:**
- `DAG.to_dict()` converts the DAG to a plain dict in the serialized form defined by `dag.md`.
- Interface metadata, step bindings, and outputs are emitted in fully explicit form so transpilers do not need to reconstruct compiler intent.
- Written as pretty-printed JSON via `dag.write(path)`.
- The JSON can be read back via `DAG.read(path)` → `DAG.from_dict()`.

---

## Public API

**File:** `api.py`

Four high-level functions:

| Function | Description |
|----------|-------------|
| `compile_workflow(path, output_path)` | Lower + force + write DAG JSON (full compilation). |
| `force_workflow(path, files)` | Lower + force, return `DAG` object. |
| `load_workflow(path, files)` | Semantic check only, return `WorkflowCheck`. |
| `transpile_dag(dag_path, target)` | Compile DAG JSON then transpile to `"cwl"`, `"wdl"`, or `"nf"`. |

---

## Pre-run-time Validation

**File:** `semantic/wf/validate.py`

Validates concrete workflow inputs against the inferred input type:
- Records must have all required fields.
- Table columns must be arrays of equal length.
- Scalars are checked for Python type (`str`, `int`, `float`).

---

## Type Utilities

**File:** `types.py`

Provides `normalize_swl_type()`, `is_optional_type()`, `is_array_type()`, `to_cwl_type()`, `to_wdl_type()`, `to_nf_qualifier()` for converting between SWL type strings and target-platform type representations.

---

## Transpilation

**Directories:** `transpile/cwl/`, `transpile/wdl/`, `transpile/nf/`

Transpilers consume the compiled DAG JSON and convert it into target-specific workflow documents. They operate on the normalized DAG contract rather than on source `.swl` or `.sh` files.

| Target | File | Format |
|--------|------|--------|
| CWL | `transpile/cwl/emit.py` | CWL v1.0 packed (`$graph`). Handles sub-workflows, record bindings (via `ExpressionTool`), `map_by` (via JS-based grouping tools), scatter, resource requirements. |
| WDL | `transpile/wdl/emit.py` | WDL v1.1. Handles struct collection, sub-workflows, scatter for mapped steps, `collect_by_key` for `map_by`, resource requirements. |
| Nextflow | `transpile/nf/emit.py` | Nextflow DSL2. Handles processes with directives, mapped steps via channels, `map_by` via `groupTuple`. |

Shared utilities live in `transpile/common.py` (identifier normalization, field chains, interpolation, run values, table column helpers).

---

## Debug / Evaluation Scripts

**Directory:** `eval/`

Six standalone CLI tools, each runnable as `python -m swl.eval.<name>`:

| Script | Purpose |
|--------|---------|
| `syntax_wf.py` | Show workflow tokens and AST tree. |
| `syntax_task.py` | Show task annotation parse result. |
| `semantic_wf.py` | Show semantic workflow view (imports, errors, inferred inputs, signature). |
| `semantic_task.py` | Show task semantic signature. |
| `ir.py` | Lower a `.swl` file to IR and print the tree. |
| `dag.py` | Force a `.swl` file to DAG and print as JSON. |

---

## Relation to the Spec

| Spec section | Implementation |
|---|---|
| Pre-condition (whitespace normalization, indent tokens) | `syntax/wf/lexer.py` |
| Workflow grammar (BNF) | `syntax/wf/parser.py` |
| Comments | `syntax/wf/lexer.py` (skipped during lexing) |
| Block syntax | `syntax/wf/parser.py:_parse_block`, `_parse_block(inner=True)` |
| Precedence | `syntax/wf/parser.py` method ordering |
| Task annotations | `syntax/task/parser.py`, `semantic/task/type.py` |
| String interpolation | `syntax/task/interpolation.py` |
| Import verification | `semantic/wf/imports.py:_load_imports` |
| Name binding rules | `semantic/wf/scope.py:check_scope` |
| Type compatibility (chaining) | `semantic/task/type.py:TypeChecker.check_chain` |
| Pipeline desugaring | `ir/lower.py:_chain_to_lambda_block` |
| Record update semantics | `ir/lower.py:_normalize_update`, `dag/merge.py` |
| map / map_by semantics | `syntax/wf/builtins.py`, `dag/evaluator.py:_force_map` |
| Workflow well-formedness | `semantic/wf/infer.py:_eval_expr` |
| Compile-time checks (1-4) | `semantic/wf/check.py:Checker.load` |
| Pre-run-time checks (5) | `semantic/wf/validate.py` |
| DAG graph | `dag/node.py:DAG`, `dag/evaluator.py` |

---

## Source files

### Core compiler
- `api.py` — Public API (compile, force, load, transpile)
- `compile.py` — CLI entry point
- `loader.py` — File cache and parsing memoization
- `repl.py` — Interactive REPL
- `types.py` — Type normalization and cross-platform type conversion

### Workflow syntax
- `syntax/wf/lexer.py` — Workflow lexer
- `syntax/wf/parser.py` — Workflow parser
- `syntax/wf/node.py` — Workflow AST node types
- `syntax/wf/builtins.py` — Shared pattern-matching for `import`, `map`, `map_by`

### Task syntax
- `syntax/task/parser.py` — Task annotation parser
- `syntax/task/bash.py` — Bash script parser
- `syntax/task/interpolation.py` — String interpolation parser

### Workflow semantic analysis
- `semantic/wf/check.py` — Workflow semantic checker (orchestrator)
- `semantic/wf/type.py` — Workflow type system (ScalarType, RecordType, TableType, FunctionType)
- `semantic/wf/imports.py` — Import resolution (.sh and .swl)
- `semantic/wf/scope.py` — Scope checking and chain type checking
- `semantic/wf/infer.py` — Symbolic evaluation for input inference
- `semantic/wf/signature.py` — Workflow signature construction
- `semantic/wf/bashvars.py` — Bash variable resolution checking
- `semantic/wf/validate.py` — Pre-run-time input validation

### Task type system
- `semantic/task/type.py` — Task type system (TaskSignature, Param, TypeKind, TypeChecker)

### Intermediate representation
- `ir/node.py` — IR node types (Literal, Input, Name, Ref, Record, Field, Update, Function, Lambda, Closure, Apply, Map, Variable, Block, Unknown)
- `ir/lower.py` — AST-to-IR lowering (chain desugaring, map/map_by detection, normalization)

### DAG construction and finalization
- `dag/node.py` — DAG data types and JSON serialization (DAG, StepCall, Input, Literal, Record, TableSource, Field, Merge, Output, OutputSpec, ForcedFunction)
- `dag/context.py` — ForceEnv variable scoping
- `dag/binding.py` — DAG binding serialization/deserialization
- `dag/evaluator.py` — IR forcing engine (force_value, _force_apply, _force_map)
- `dag/forcer.py` — Top-level forcing orchestration (ForceState, force, force_file)
- `dag/emit.py` — Step call emission (task, workflow, mapped step)
- `dag/merge.py` — Merge canonicalization and flattening
- `dag/finalize.py` — DAG finalization (metadata, pruning, output specs)
- `dag/tooldefs.py` — Tool definition construction and workflow materialization

### Transpilers
- `transpile/common.py` — Shared transpiler utilities
- `transpile/_cli.py` — Shared CLI argument parsing
- `transpile/cwl/emit.py` — CWL v1.0 transpiler
- `transpile/wdl/emit.py` — WDL v1.1 transpiler
- `transpile/nf/emit.py` — Nextflow DSL2 transpiler

### Debug/evaluation
- `eval/syntax_wf.py` — Show workflow tokens and AST
- `eval/syntax_task.py` — Show task annotation parse
- `eval/semantic_wf.py` — Show semantic workflow view
- `eval/semantic_task.py` — Show task semantic signature
- `eval/ir.py` — Show lowered IR tree
- `eval/dag.py` — Show forced DAG
