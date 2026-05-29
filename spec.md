# SWL Specification

## Principles

- Minimalistic design: Specialized for building workflows
- Incremental build: Users can build workflows piece by piece
- Fail fast: Catch as many issues as possible during static time
- Efficient execution: Lazy evaluation and caching during run time

## Syntax

### GBNF Grammar (Workflow Language)

#### Pre-condition
- Any "\r\n" is replaced with "\n". Any "\r" is removed.
- "\n" is preserved.
- each block of indents (4 spaces) has been replaced with `bstart` and `bend`
  tokens that represent start (indent) and end (de-indent) of an indented block
- all other non-newline whitespace has been removed

```bnf
root           ::= block (eol block)*

block          ::= (binding eol)* expr

binding        ::= name "=" expr

expr           ::= value | operation | lambda

lambda         ::= "\" name "->" block

operation      ::= apply | parens | get | update | chain

parens         ::= "(" expr ")"

apply          ::= name ws expr

get            ::= name ("." name)+

update         ::= expr ( "//" expr)+

chain          ::= expr ( "|" expr)+

record         ::= "{" pairs ","? ws? "}"

value          ::= name | number | string | record

pairs          ::= pair?  |  pair (ws? "," ws? pair)+

pair           ::= name ws? ":" ws? expr

name           ::= [a-zA-Z_][a-zA-Z0-9_]*

number         ::= [0-9]+  "."?  [0-9]*

string         ::= "\"" ([^"\] | "\\" . )* "\""

eol             ::= "\n"

ws             ::= [ \t\n]+
```

#### Comments

Comments are prefixed by `#`

#### Block Syntax (described in words)

A block is one or more statements, each ending with `\n` (newline).

- A **lambda body** is either:
  - An **indented block**: multiple lines where all lines after the first are indented more than the `->` line
  - A **single expr**: on the same line as `->`

- An **indented block** ends when:
  - A line is encountered with indentation less than or equal to the `->` line, OR
  - End of file is reached

- The **final line in a block must be an `expr`** (not a binding).

#### Types

- **record** (`rec`): A collection of key-value pairs where values are scalars, functions, or other records: `{ str: str|file|num|rec, ... }`
- **table** (`tab`): A specialized record whose top-level fields are homogeneous arrays of equal length: `{ str: [t], ... }`
- A `tab` is the canonical batch value in SWL. It is record-shaped, but type-distinct from plain `rec` for workflow typing.
- If a table contains fields with uneven array lengths, it triggers a static validation error.
- Row semantics for a table are derived by zipping aligned array positions across all fields; rows are not a separate surface type.

- **task**: `rec -> rec`
- **simple workflow**: `rec -> rec`
- **batch workflow**: `tab -> tab` or `tab -> rec`

#### Built-in functions

- `import: str -> fun` imports a task or workflow as a function
- `map: (rec -> rec) -> tab -> tab` maps a `rec -> rec` function over the logical rows of a table and returns a table
- `map_by: (rec -> rec) -> str -> tab -> tab` maps a function over table, grouping the rows
  of the table by a key, to return a table

#### Workflow well-formedness

- A workflow must evaluate to a function.
- A workflow whose final value is a scalar, record, or saturated computation is invalid.

#### Name binding rules

- Variable names are unique within their immediate scope.
- A variable name cannot be bound twice in the same scope.
- Import bindings follow the same rule; duplicate import names in the same scope are errors.
- Any nested scope may shadow names from outer scopes.
- This shadowing rule applies uniformly to lambda parameters and bindings inside nested blocks.


#### Precedence (highest to lowest)

1. Record get (`.`)
2. Function application (`ws`)
3. Record update (`//`)
4. Lambda arrow (`->`)
5. Function chain (`|`)
6. Let binding (`=`)

---

### GBNF Grammar (Annotation Language)

Annotations are comment block in bash scripts, starting with `#`.

#### Assumptions
- The prefix `#` with an optional " " is removed.
- Any "\r\n" is replaced with "\n". Any "\r" is removed.
- "\n" is preserved

```bnf
annotation    ::= task-doc section+

task-doc      ::= "@" phrase eol+

section       ::= "in"  ws param-list
                | "out" ws param-list
                | "run" ws param-list

param-block   ::= param | param param-block

param         ::= names sp? type? sp? default? sp? desc* eol+

names         ::= name ("," sp? name)*

type          ::= (simple-type | array-type) "?"?

simple-type   ::= "file" | "str" | "int" | "float"

array-type    ::= "[file]" | "[str]" | "[int]" | "[float]"

default       ::= "=" value

value         ::= word | string

#### Task parameter scope rules
- Input parameters form one scope.
- Output parameters form one scope.
- Run parameters form one scope.
- Duplicate names within a single section scope are errors.
- It is legal for an input parameter to have the same name as an output parameter.

desc          ::= eol? "|" phrase

phrase        ::= [^\n]+

word          ::= [^\n\t |]+

string        ::= "\"" [^"\n]* "\""

name          ::= [a-zA-Z_][a-zA-Z0-9_]*

eol            ::= "\n"

sp            ::= [ \t]+

ws            ::= [ \t\n]+
```

#### Features

- Multiple names per line: `fastq1, fastq2 file`
- Default value binding: `= 2`
- String interpolation: `${outbase}.bam`
- Descriptions: prefixed by `|` can occur on the same line or on continuation lines

---

## Type System

### Types

| Type    | Description     |
|---------|-----------------|
| `file`  | File path       |
| `str`   | String          |
| `int`   | Integer         |
| `float` | Floating point  |

### Type Annotations (for documentation)

In task scripts:

```
# in  fastq1 file                  | description
# out bam    file = ${outbase}.bam | description
# run cpu         = 2
```

- `# in`  - input parameter (required)
- `# out` - output parameter (required)
- `# run` - runtime resource
- `?` postfix marks optional input parameter: `file?`

---

## Compile-time Checks

Upon compilation of the workflow ...

### 1. Import Verification

- Imported files must exist at compile time
- Import paths are relative to the location of the workflow file
- Circular imports are a compile error

### 2. Type Compatibility

**Pipeline (`a | b`):**
- For each output of `a` with the same name as an input of `b`,
  the types must match, and the input of `b` may be optional

**Record Union (`r1 // r2`):**
- Same field name: types must match exactly
- Right record overrides left (feature)

**Type Matrix:**

If functions `a` and `b` are chained together,
any record fields with the same name must have compatible types:

| a.out  | b.in     | Allowed? |
|--------|----------|---------|
| `file` | `file`   | ✓       |
| `file` | `file?`  | ✓       |
| `int`  | `int`    | ✓       |
| `int`  | `int?`   | ✓       |
| `str`  | `str`    | ✓       |
| `str`  | `str?`   | ✓       |

All other combinations are not allowed.

### 3. Workflow function and DAG

- Workflow must return a function
- DAG graph will be constructed and checked for circularity

### 4. Inference of Workflow Input

- Workflow input type is inferred from the unsatisfied inputs of all tasks
- A task input is satisfied if there is an upstream task output with the same
  name


## Pre-run-time Checks

Upon providing inputs to a workflow ...

### 5. Workflow Input Requirements

- Provided inputs will be type-checked against the workflow inputs
- Every required input (`type` without `?`) must be satisfied by
  an workflow input or be satisfiable by an output from an upstream task


## Run-time Checks

### 6. Dependency

- Return code of each task will be checked
- Downstream task will only be run if all of its dependent tasks succeed


## Post-run-time Checks

### 7. Completion Status

- Return code of final task
- Existence of each workflow output

---

## Semantics

### Records
- Collection of key-value pairs
- Keys are field names
- Values are numbers, strings, records, functions, or arrays where the type system permits them

### Tables
- A table is the canonical batch value
- A table is columnar: each field is an array, and all top-level arrays have the same length
- A table has derived row semantics: row `i` is formed by taking element `i` from each field array
- Field access on a table is ordinary record field access; if `xs : tab` and `field : [t]`, then `xs.field : [t]`
- SWL does not define a separate array-of-record batch type for workflow semantics

### Functions
- Importing a task or a workflow returns a function
- Functions can defined using a lambda

### Function Application
- `f r` applies function `f` to a record or table argument of the required type
- Partial application returns a new function

### Map
- `map` is a builtin
- If `f : rec -> rec`, then `map f : tab -> tab`
- `map f xs` applies `f` to each logical row of `xs`
- The result is reassembled as a table
- If `f` is batch-typed, then `map f` is a compile-time error for now

### Pipeline
- `A | B | C` desugars to:
```
\x ->
    a = A x
    b = B (x // a)
    c = C (x // a // b)
    a // b // c
```

### Record Update

- `r1 // r2` merges records or tables, right overrides left
- Extra fields are allowed (ignored)
- If a table `t` is updated with a record `r` via `t // r`, scalar properties in `r` are implicitly duplicated across all rows of `t` to preserve table length integrity. Similarly for `r // t`.
- If both sides are tables, matching fields must remain array-typed and length-compatible.

---

## Examples

See `tests/` directory.

