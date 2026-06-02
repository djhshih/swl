# SWL Architecture

SWL is a workflow language that compiles to a logical DAG of task invocations. The pipeline has these stages:

1. **Workflow Lexing** → token stream
2. **Workflow Parsing** → AST (`syntax/wf/node.py`)
3. **Task Parsing** → parsed task annotation + bash body (`syntax/task/node.py`)
4. **Semantic Checking** → validated imports, scope, types, inferred signature (`semantic/wf/check.py`)
5. **IR Lowering** → semantic IR tree (`ir/node.py`)
6. **IR Forcing** → concrete DAG (`ir/dag.py`)
7. **Output** → JSON-serialized DAG emission

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
| `chain`  | `Chain`      | `left \| right` (pipeline) |

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

**File:** `semantic/wf/check.py`

**Input:** Parsed workflow AST (from stage 2) plus imported tasks/workflows.

**Output:** A `WorkflowCheck` object containing: validated AST, imports, errors, inferred inputs, workflow signature, workflow type, batch-vs-simple classification.

**What happens:**

1. **Import resolution** (`_load_imports`): Scans top-level bindings for `name = import("path")` calls. `.sh` paths → task imports; `.swl` paths → workflow imports (recursively loaded). Circular imports are detected via a `_loading` stack.

2. **Scope checking** (`_check_scope`): Walks the AST ensuring no duplicate variable names within the same scope. Nested scopes (lambda params, inner blocks) may shadow outer names.

3. **Chain type checking** (`_check_chains`): Detects `A | B` chains in the AST and calls the task-level `TypeChecker.check_chain()` for any pair of known tasks.

4. **Input inference** (`_infer_inputs`): Evaluates the workflow body in a symbolic/semantic fashion to determine what inputs the workflow demands. Uses a set of semantic value types:
   - `OpenRecord` / `ClosedRecord` — record values with known/unknown fields.
   - `FunctionValue` — a known function with a signature.
   - `ClosureValue` — a partially applied function.
   - `ComputationValue` — a fully applied function (a task invocation).
   - `TableValue` — a table (columnar batch data).
   - `UnknownValue` — unresolved.
   - `TypedValue` — a value with a known type but unknown origin.

   The semantic evaluator (`_eval_expr`) recursively walks the AST. When it encounters a field access like `x.field` on an `OpenRecord`, it adds `field` to the `demanded` set (this is how input fields are inferred). When it encounters a function application, it checks saturation — if all required inputs are available, it produces a `ComputationValue`; otherwise a `ClosureValue`.

5. **Signature construction** (`_build_workflow_signature`): Determines whether the workflow is batch (`tab -> tab`) or simple (`rec -> rec`) based on whether `map`/`map_by` is used. Builds the final `TaskSignature` from inferred inputs and outputs. Also produces `workflow_type` as a `FunctionType(input_type, output_type)` from `semantic/wf/type.py`.

**Type AST** (`semantic/wf/type.py`):
- `ScalarType("file" | "str" | "int" | "float" | "?")`
- `ArrayType(item)`
- `RecordType(fields, open)`
- `TableType(columns)`
- `FunctionType(input, output)`

---

## Stage 6: IR Lowering

**Files:** `ir/lower.py`, `ir/node.py`

**Input:** Workflow AST + imports + signature from the semantic checker.

**Output:** A semantic IR tree of `ir.Node` dataclass instances.

**IR node types:**
| Node       | Meaning |
|------------|---------|
| `Literal`  | A literal value (string, number) |
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
- **Lambda s** become `ir.Lambda` with a `Block` body.
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

## Stage 7: IR Forcing

**Files:** `ir/force.py`, `ir/dag.py`

**Input:** Semantic IR tree (from stage 6).

**Output:** A concrete `DAG` object with resolved inputs, steps, and outputs, ready for JSON serialization.

**Key DAG types** (`ir/dag.py`):
| Type | Meaning |
|------|---------|
| `Input` | A named workflow input with optional type/desc/optionality |
| `Literal` | A literal value |
| `Record` | A record of named DAG values |
| `TableSource` | A table source (input columns) |
| `Field` | Projection of a field from a step or input |
| `Merge` | A transient merge of two DAG values (left // right) during forcing |
| `StepCall` | A concrete task or workflow invocation |
| `MappedStep` | A task invocation wrapped in map/map_by |
| `Output` | A named workflow output |
| `ForcedFunction` | A partially-saturated function reference |
| `DAG` | The top-level container: `inputs`, `steps`, `outputs` |

**What happens:**

The `Forcer` "forces" the IR by resolving all function applications into concrete step calls:

1. **`force_value(node, env)`** — Recursively evaluates IR nodes in a `ForceEnv`:
   - `Literal` → `dag.Literal`
   - `Name` → either a `ForcedFunction` (if it's a known variable bound to a function) or an `Input` (workflow input)
   - `Record` → `dag.Record` with forced fields
   - `Field` → `dag.Field` with attempted projection
   - `Update` → `dag.Merge` or record merge
   - `Function` / `Lambda` → `ForcedFunction` wrapper
   - `Apply` → tries to apply the function to the argument
   - `Map` → maps function over table source
   - `Block` → evaluates in a local `ForceEnv`
   - `Variable` → forces the value and caches it

2. **Function application** (`_apply`):
   - If the function is a built-in (`map`, `map_by`), handles partial application (binding `f`, then `key`/`xs`).
   - If the argument is a `MappedStep`, delegates to `_apply_mapped`.
   - For tasks: checks saturation — if all required inputs are available in the bound record, emits a `StepCall` via `_emit_task_call`. Otherwise returns a `ForcedFunction` (partial application).
   - For workflows: similar saturation check, emits `StepCall` with `type='workflow'`.
   - For lambdas: if the bound argument is a `Record`, evaluates the lambda body directly.

3. **Step call emission** (`_emit_task_call` / `_emit_workflow_call`):
   - Normalizes task inputs by projecting named fields from the bound argument.
   - Creates a `StepCall` or `MappedStep` with unique ID, dependencies, and tool definition.
   - Caches step calls by function path + input keys to deduplicate.

4. **Map emission** (`_emit_mapped_step`):
   - Creates a `MappedStep` with map source information (table columns, grouping key for `map_by`).
   - Sets `input_schema` / `output_schema` from the function signature.

5. **Root forcing** (`_force_root`):
   - If the final value is a `ForcedFunction`, applies it to synthetic input placeholders to saturate it.
   - For batch workflows, creates a table input placeholder and applies `map`.

6. **DAG finalization** (`_finalize_dag`):
   - `_refine_input_metadata()`: Materializes workflow input metadata by merging inferred interface information with task input specs, including type, description, and optionality where available.
   - `_final_outputs()`: Converts the forced root value into the explicit top-level workflow outputs required by `dag.md`, collecting top-level record fields as named outputs or wrapping a single value as `result`.
   - `_prune_unused_inputs()`: Removes input ports that aren't actually referenced by any step or output.
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

**File:** `compile.py` (entry point), `ir/dag.py` (serialization)

**Input:** A normalized `DAG` object.

**Output:** A JSON file written to `dag/<name>.json`.

**What happens:**
- `DAG.to_dict()` converts the DAG to a plain dict in the serialized form defined by `dag.md`.
- Interface metadata, step bindings, and outputs are emitted in fully explicit form so transpilers do not need to reconstruct compiler intent.
- Written as pretty-printed JSON via `dag.write(path)`.
- The JSON can be read back via `DAG.read(path)` → `DAG.from_dict()`.

---

## Pre-run-time Validation

**File:** `semantic/wf/validate.py`

Validates concrete workflow inputs against the inferred input type:
- Records must have all required fields.
- Table columns must be arrays of equal length.
- Scalars are checked for Python type (`str`, `int`, `float`).

---

## Transpilation

**Directories:** `transpile/cwl/`, `transpile/wdl/`, `transpile/nf/`

Transpilers consume the compiled DAG JSON and convert it into target-specific workflow documents. They operate on the normalized DAG contract rather than on source `.swl` or `.sh` files.

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
| Import verification | `semantic/wf/check.py:_load_imports` |
| Name binding rules | `semantic/wf/check.py:_check_scope` |
| Type compatibility (chaining) | `semantic/task/type.py:TypeChecker.check_chain` |
| Pipeline desugaring | `ir/lower.py:_chain_to_lambda_block` |
| Record update semantics | `ir/lower.py:_normalize_update`, `ir/force.py:_merge_values` |
| map / map_by semantics | `syntax/wf/builtins.py`, `ir/force.py:_force_map` |
| Workflow well-formedness | `semantic/wf/check.py:_infer_inputs` |
| Compile-time checks (1-4) | `semantic/wf/check.py:Checker.load` |
| Pre-run-time checks (5) | `semantic/wf/validate.py` |
| DAG graph | `ir/dag.py:DAG`, `ir/force.py` |


## Source files

- syntax/wf/lexer.py — Workflow lexer
- syntax/wf/parser.py — Workflow parser
- syntax/wf/node.py — Workflow AST node types
- syntax/wf/builtins.py — Shared pattern-matching for `import`, `map`, `map_by` (used by checker and lowerer)
- syntax/task/parser.py — Task annotation parser
- syntax/task/bash.py — Bash script parser
- syntax/task/interpolation.py — String interpolation parser
- semantic/wf/check.py — Workflow semantic checker (imports, scope, type inference, partial application)
- semantic/wf/type.py — Workflow type system (ScalarType, ArrayType, RecordType, TableType, FunctionType)
- semantic/wf/validate.py — Pre-run-time input validation
- semantic/task/type.py — Task type system (TaskSignature, Param, TypeKind, type compatibility matrix)
- ir/node.py — IR node types (Function, Lambda, Closure, Apply, Map, Block, Variable, etc.)
- ir/lower.py — AST-to-IR lowering (chain desugaring, map/match_map_by, normalization)
- ir/force.py — IR forcing (DAG construction, step call emission, mapped step emission, partial application, root forcing)
- ir/dag.py — DAG dataclasses and JSON serialization (StepCall, MappedStep, Input, Literal, Field, Merge, Record, TableSource, DAG)
- transpile/cwl/emit.py — CWL transpiler (reference pattern for new transpilers)
- transpile/cwl/test_emit.py — CWL transpiler tests
- compile.py — Compiler entry point
- repl.py, eval\*.py — Debug evaluation scripts

