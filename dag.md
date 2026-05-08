# DAG-oriented semantic IR redesign

## Goal
Make workflow lowering produce a scope-correct, DAG-friendly semantic IR so `force.py` can evaluate semantic identities instead of reconstructing sharing from copied trees.

## Core canonical shape
Every workflow lambda should lower to:

```python
Lambda(param, Block(bindings=[Variable(...), ...], result=...))
```

This is mandatory even when there are no intermediate bindings.
An empty body is represented as:

```python
Lambda(param, Block(bindings=[], result=...))
```

## Semantic nodes

### `Variable`
A bound semantic computation.

Suggested shape:
- `id: int`
- `name: str`
- `value: Node`

A source binding like:

```swl
a = align x
```

lowers to one `Variable`.

### `Ref`
A use of a previously defined semantic variable.

Suggested shape:
- `id: int`
- `name: str`

A later use of `a` should become `Ref(id_of_a, 'a')`, not a copied `Apply(...)` tree.

### `Block`
A let-graph body:
- `bindings: List[Variable]`
- `result: Node`

### `Lambda`
A workflow function:
- `param: str`
- `body: Block`
- `signature: Optional[TaskSignature]`

Every lambda body should be a `Block`.

## Lowering rules

### Local names
When a name resolves to a local workflow binding:
- return `Ref(id, name)`
- do not inline the bound value

### Imports
Imported tasks/workflows remain `Function` values in the lowering env and are not emitted as `Variable`s.

### Bindings
When lowering a source binding:
1. allocate a fresh variable id
2. lower the RHS once
3. emit `Variable(id, name, rhs)`
4. update the env so later references to that name become `Ref(id, name)`

### Lambda bodies
Every lambda body should be wrapped in a `Block`, even if the source body has no bindings.

## Chain normalization
`Chain` should not survive as a canonical post-normalization form.
Instead, chain syntax should normalize to the same explicit staged block graph as function syntax.

Example:

```swl
align | sort | call
```

should normalize approximately to:

```python
Lambda(
  param='_input',
  body=Block(
    bindings=[
      Variable(v1, '_s1', Apply(align, Name('_input'))),
      Variable(v2, '_s2', Apply(sort, Update(Name('_input'), Ref(v1, '_s1')))),
      Variable(v3, '_s3', Apply(call, Update(Update(Name('_input'), Ref(v1, '_s1')), Ref(v2, '_s2')))),
    ],
    result=Update(Update(Ref(v1, '_s1'), Ref(v2, '_s2')), Ref(v3, '_s3')),
  ),
)
```

So chain sugar and explicit staged workflows converge to the same semantic representation.

## Force-time model
`force.py` should evaluate `Block` + `Variable` + `Ref` as a semantic graph.

### Main rule
Each `Variable.id` is forced at most once.

### Forcing behavior
- register all variables visible in the current block
- when a `Ref(id, ...)` is encountered:
  - if already forced, reuse it
  - otherwise force the corresponding variable once and cache it

This makes DAG sharing explicit and scope-correct.

## Why this is better
- no duplicated semantic work
- scope is encoded by unique variable ids
- `function.swl` and `pipe.swl` converge before forcing
- `force.py` becomes a let-graph evaluator rather than a tree-reconstruction pass
