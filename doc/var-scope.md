# Variable Scope

## Shared `Scope` class

All compiler phases share a single `Scope` class defined in `python/swl/semantic/scope.py`:

```python
class Scope:
    def __init__(self, parent=None):
        self.parent = parent
        self.locals = {}               # name -> Binding

    def declare(self, name, ...):
        if name in self.locals:
            raise DuplicateBindingError(name)  # same-scope only
        self.locals[name] = Binding(...)

    def set_local(self, name, ...):
        # set field on existing binding, or declare new one
        ...

    def resolve(self, name):
        if name in self.locals:
            return self.locals[name]
        if self.parent is not None:
            return self.parent.resolve(name)
        return None                     # undefined
```

Key design: `declare()` checks only `self.locals` (same-scope duplicates). Names in parent scopes are NOT checked, so shadowing across scopes is allowed. `resolve()` walks the parent chain for lookup.

Each `Block` and `Lambda` creates a child `Scope(parent=current)`. Blocks are NOT flat copies of the parent set — they inherit via the parent pointer, so adding a name to a child scope never affects parent scopes.

---

## Annotation Language (Task Scripts)

### Scope regions

All three annotation sections (`in`, `out`, `run`) share a single flat scope. 

**Cross-section duplicate rules** (implemented in `signature_from_task` at `type.py:170`):

| Sections | in+out | in+run | out+run |
|----------|--------|--------|---------|
| Shadowing allowed? | Yes | No | No |

- `in x` + `out x` is **allowed** (common pattern: sort takes `bam`, produces `bam`).
- `in x` + `run x = ...` is an error.
- `out x` + `run x = ...` is an error.
- Same-section duplicates (`in x` + `in x`, etc.) are always errors.

The implementation uses three separate `set` accumulators (`in_names`, `out_names`, `run_names`). The `in` section checks against `in_names` (duplicate) and `run_names` (cross-section). The `out` section checks against `out_names` and `run_names`. The `run` section checks against all three.

### Interpolation scope

Interpolations (`${name}`, `$name`) resolve against the following sets:

| Location | Resolution scope |
|----------|-----------------|
| Output default (`= ${outbase}.bam`) | Input parameters + run parameters |
| Run parameter default (`= ${some_var}`) | All parameters (the unified annotation scope) |
| Command block body | Input parameters + run parameters + preceding bash assignment LHS + shell builtins |

Bash builtins are recognized via an allowlist in `python/swl/syntax/task/bash.py` (`_BUILTIN_VARS`) and excluded from validation.

The variable checker in `python/swl/semantic/wf/bashvars.py` walks bash `Script` AST statements sequentially, maintaining a mutable `defined` set. On each `Assignment`, the LHS is added to `defined`. On each `Command`, all RHS variable references are checked against `defined`.

---

## Workflow Language

### Grammar structure

A lambda's grammar is `\param -> block`. The lambda introduces the parameter into scope via a child `Scope`, then its body is a block. The resulting scope chain:

```
root Scope → outer block Scope → lambda Scope → inner block Scope
```

A block (`(binding eol)* expr`) creates a child `Scope`. Bindings add names to that block's scope. The final line must be an expression. Blocks nest:

```
\x ->                  # lambda Scope: {x}
    y = 1              #   inner block Scope: {y}; parent has {x}
    \z ->              #   lambda Scope: {z}; parent has {x, y}
        y = 2          #     innermost block Scope: {y}; shadows outer y
        y + z          #     final expr
```

### Shadowing rules

All four compiler phases now use the same `Scope` class with a consistent policy:

| Phase | Shadowing policy | Mechanism |
|-------|-----------------|-----------|
| Scope checker (`walk_scope`) | **Allowed** across scopes, **error** same-scope | `Scope.declare()` checks only `self.locals`. Block creates `Scope(parent=...)`. Duplicate in same block → error. Shadowing across blocks → allowed. |
| Inference (`infer.py`) | **Allowed** | `Scope.set_local()` overwrites the binding value. |
| IR lowering (`lower.py`) | **Allowed** | `Scope.set_local()` stores the IR node. |
| DAG forcing (`evaluator.py`) | **Allowed** | `Scope.resolve()` returns the innermost match. |

### Forward reference detection

Before processing binding `i` in a block, the scope checker builds `defined_so_far` by flattening all ancestor scopes' names plus names of bindings 0..i-1. It collects `Id` references in binding `i`'s value expression and checks:

1. Ref equals binding's own name → self-reference error.
2. Ref not in `defined_so_far` but matches a later binding's name → forward reference error.
3. Ref not in `defined_so_far` and not a later binding → free variable (may resolve later as workflow input).

`collect_name_refs` skips lambda bodies (`NodeType.fun`) to avoid false-positive forward references inside nested functions.

### Import bindings

An import binding (`name = import "path"`) is syntactically a regular binding. At the top level, imports are pre-resolved by `load_imports` and stored in the workflow's `imports` dict. Inside a lambda body, `import "path"` is handled inline via `match_import` and `_lower_inline_import`.

The `Scope.resolve()` call in `lower_binding` has a guard that prevents the top-level-import short-circuit from firing on shadowed names: if the name is already in scope (via the parent chain), it falls through to normal expression evaluation.

---

## Current architecture per compiler phase

### Phase 1: Scope checker (`python/swl/semantic/wf/scope.py`)

**Data structure**: `Scope` class with parent pointer.

**Two-pass design**:

**Pass 1 — `walk_scope`**: Walks the AST to build the `Scope` tree and detect errors:
- For a **block**: creates `Scope(parent=outer_scope)`, iterates body items. For each binding, calls `scope.declare(name)` (catches same-scope duplicates). Binding values are walked with the same scope (so they see all prior bindings). Import bindings skip the value walk (`match_import` check).
- For a **lambda/fun**: creates `Scope(parent=outer_scope)`, declares the param name, recurses into the body.
- After registering bindings, runs forward-reference detection on each binding's value expression.

**Pass 2 — `annotate_ast`** (`_walk_annotate`): After `walk_scope` builds the tree and all names are registered, this pass walks the AST again and annotates every `Id` node with its resolved `Binding`:

```python
def _walk_annotate(checker, expr, scope):
    if expr.type == wf_node.NodeType.id:
        binding = scope.resolve(expr.name)
        if binding is not None:
            expr.binding = binding        # <-- Id.binding populated here
        return
    # recurse into children with appropriate Scope
```

The `Identifier.binding` field is defined in `python/swl/syntax/wf/node.py`.

### Phase 2: Input inference (`python/swl/semantic/wf/infer.py`)

**Data structure**: `Scope` frames storing inferred semantic values (`FunctionValue`, `OpenRecord`, `ClosureValue`, etc.) on `Binding.value`.

**Mechanism**:
- For a **block**: creates `Scope(parent=current)`. Bindings are evaluated in order; results stored via `scope.set_local(name, value=result)`.
- For a **lambda**: creates `Scope(parent=current)`, binds the parameter name. The resulting `FunctionValue` captures the current scope for closure semantics.
- For an **`id`** node: reads `id_node.binding.value` (populated by the annotation pass) — no dict lookup needed.

Import bindings are skipped in `eval_prefix_bindings` (they are pre-resolved in the `imports` dict). This function is only called for top-level prefix bindings.

### Phase 3: IR lowering (`python/swl/ir/lower.py`)

**Data structure**: `Scope` frames storing IR nodes (`ir.Ref`, `ir.Name`, `ir.Input`, `ir.Function`) on `Binding.ir_node`.

**Mechanism**:
- For a **block**: creates `Scope(parent=current)`. For each binding: lower the value expression, create `ir.Variable(id, name, value)`, store `ir.Ref(id, name)` via `scope.set_local(name, ir_node=...)`.
- For a **lambda**: creates `Scope(parent=current)`, declares param name with `ir_node=ir.Name(name)`. Body is lowered with this scope.
- For an **`id`** node: `_lower_name` calls `scope.resolve(name)` and reads `binding.ir_node`. Falls through to `imports[name]` or `signature.inputs[name]` for names not in scope.
- **`lower_binding`**: has a short-circuit for top-level imports guarded by `scope.resolve(name) is None`. This ensures the short-circuit only fires for first-use top-level imports, not shadowed names in inner scopes.

**Closure captures**: Computed by `_compute_lambda_captures`, which scans the lambda body for free variable references to outer scope entries. Captures are stored on `ir.Lambda.captures` as `{name: ir_node}` dict.

### Phase 4: DAG forcing (`python/swl/dag/evaluator.py` + `python/swl/dag/forcer.py` + `python/swl/dag/tooldefs.py`)

**Data structure**: `Scope` frames replacing the old `ForceEnv` chain. A parallel `forcer.variables` dict maps `int` (variable ID) → `ir.Variable` for ID-based resolution of `ir.Ref` nodes.

**Mechanism**:
- For a **block**: creates `Scope(parent=env)`. Pre-registers all bindings in `forcer.variables` by ID. Binds each name to `ir.Ref(id, name)` in the new `Scope`.
- For a **lambda application**: creates a `Scope()` for the lambda body. Restores captured variables from `node.captures` into the scope. Lambda params are `ir.Name` nodes resolved via `Scope.resolve()`.

**Two resolution mechanisms**:
- `ir.Name` — resolved via `Scope.resolve(name)` (lexical chain). Used for lambda parameters.
- `ir.Ref` — resolved via `forcer.variables[id]` (ID-based, bypasses chain). Used for local variable bindings.

---

## Effects of `import` inside a lambda

### Calling `import` as an expression

`import "path"` is recognized by `match_import`. It can appear anywhere an expression is expected:

- Top-level binding: `align = import "align.sh"`
- Lambda body binding: `\x -> t = import "task.sh"; t x`
- Direct application: `\x -> import "task.sh" x`

All cases detect the `apply(Identifier("import"), String(path))` pattern and load the file.

### Binding `import` result to a variable within a lambda

```
\x ->
    t = import "task.sh"
    t x
```

**Scope-checker path** (`scope.py`):
1. `walk_scope` for lambda creates `Scope(parent=outer)`, declares `x`.
2. `walk_scope` for body block creates child `Scope`, declares `t`.
3. Import value is skipped (no name refs to check).
4. `annotate_ast` annotates `t` in the body with the binding for the inner `t`.

**Inference path** (`infer.py`):
1. Lambda creates child `Scope`, binds `x`.
2. `_eval_block_items` processes binding: evaluates `import "task.sh"` → `FunctionValue`.
3. `t` bound in scope via `set_local`.

**IR-lowering path** (`lower.py`):
1. Lambda creates child `Scope`, declares `x` with `ir_node=Name("x")`.
2. `_lower_block_items` processes binding: `lower_binding` checks guard — `t` not in scope → short-circuit → `_function_from_import` → `Function`. Wait — actually `t` is NOT in the top-level imports dict (it's a new binding), so `"t" in imports` is False. Falls through to `lower_expr(expr.value, ...)` which matches `import "task.sh"` → `_lower_inline_import` → `Function`.
3. `Variable` created, `scope.set_local("t", ir_node=Ref(...))`.

### Shadowing

```
align = import "align.sh"     # top-level import
\x ->
    align = import "sort.sh"   # shadows outer align
    align x                    # refers to sort.sh, not align.sh
```

Lowering for the inner `align = import "sort.sh"`:
- `lower_binding(expr, local_scope, imports)`:
  - `"align" in imports` → True (top-level import).
  - `local_scope.resolve("align")` → walks parent chain → finds `align` in outer block scope → NOT None.
  - Guard `scope.resolve(name) is None` is **False** → short-circuit NOT triggered.
  - Falls through to `lower_expr(expr.value, ...)` → correctly lowers `import "sort.sh"` → creates `Variable` with Function(sort.sh).

