# SWL Specification

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

#### Block Syntax (described in words)

A block is one or more statements, each ending with `\n` (newline).

- A **lambda body** is either:
  - An **indented block**: multiple lines where all lines after the first are indented more than the `->` line
  - A **single expr**: on the same line as `->`

- An **indented block** ends when:
  - A line is encountered with indentation less than or equal to the `->` line, OR
  - End of file is reached

- The **final line in a block must be an `expr`** (not a binding).

#### Built-in functions
- `import` imports a task or workflow as a function


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
- Values are numbers, strings, or records

### Functions
- Importing a task or a workflow returns a function
- Functions can defined using a lambda

### Function Application
- `f record` applies function `f` to `record`
- Partial application returns a new function

### Pipeline
- `A | B | C` desugars to:
```
\x ->
    a = A x
    b = B (x // a)
    c = C (x // a // b)
    a // b // c
```

### Record Merge
- `r1 | r2` merges records, right overrides left
- Extra fields are allowed (ignored)

---

## Examples

See `tests/` directory.

