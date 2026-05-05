# SWL Specification

## Syntax

### GBNF Grammar (Workflow Language)

```bnf
root           ::= block  |  block "\n" root

block          ::= expr   |  binding "\n" expr

binding        ::= name "=" expr

expr           ::= value | operation | lambda

lambda         ::= "\" name "->" block

operation      ::= app | parens | get | union | pipeline

parens         ::= "(" operation ")"

app            ::= name value

get            ::= name "." name

union          ::= record ("|" record)+

pipeline       ::= name ("|>" "\n"? pipeline)+

value          ::= number | string | record

record         ::= "{" "\n"? pairs ","? "\n"? "}" | name

pairs          ::= pair?  |  pair ("," "\n"? pair)+

pair           ::= name ":" value

name           ::= [a-zA-Z_][a-zA-Z0-9_]*

number         ::= [0-9]+  "."?  [0-9]*

string         ::= "\"" ([^"]*) "\""
```

### Precedence (highest to lowest)

1. Function application (whitespace)
2. Record merge (`|`)
3. Pipeline (`|>`)
4. Lambda arrows (`->`)
5. Let bindings (`=`)

---

### BNF Grammar (Annotation Language)

Annotation comments appear in bash scripts as comment lines starting with `#`.

```bnf
annotation    ::= task_doc
                | section

task_doc       ::= "#@" STRING

section        ::= "#" "in" newline param_list
                | "#" "out" newline param_list
                | "#" "run" newline param_list

param_list     ::= param
                | param param_list

param          ::= "#" NAME_LIST TYPE default_opt DESC

NAME_LIST     ::= NAME
                | NAME "," NAME_LIST

default_opt   ::= EMPTY
                | "=" VALUE

DESC          ::= EMPTY
                | "|" LITERAL
                | newline "|" LITERAL

NAME          ::= [a-zA-Z_][a-zA-Z0-9_]*
TYPE          ::= "file" | "file" "?"
              | "str" | "str" "?"
              | "int" | "int" "?"
              | "float" | "float" "?"
VALUE         ::= ${ NAME } | LITERAL
```

| Token      | Pattern                    |
|------------|----------------------------|
| EMPTY      |                            |
| NAME       | `[a-zA-Z_][a-zA-Z0-9_]*`   |
| NEWLINE    | `\n`                       |
| LITERAL    | `[^\n]+`                   |


Note: The annotation language is embedded in bash scripts as comments.
The lexer extracts comment lines and parses the annotation directives.

---

## Type System

### Types

| Type    | Description       |
|---------|------------------|
| `file`  | File path        |
| `str`   | String          |
| `int`   | Integer         |
| `float` | Floating point  |

### Type Annotations (for documentation)

In task scripts:

```
# in  fastq1   file    | description
# out bam      file   = ${outbase}.bam | description
# run cpu      int    = 2
```

- `# in`  - input parameter
- `# out` - output parameter
- `# run` - runtime resource
- `?` postfix marks optional input parameter: `file?`

---

## Compile-Time Checks

### 1. Import Verification
- Imported file must exist at compile time
- Circular imports are a compile error

### 2. Required Inputs
- Every required input (`type` without `?`) must be satisfied by:
  - Workflow input (`argv`)
  - Output from upstream task in pipeline

### 3. Type Compatibility

**Pipeline (`a |> b`):**
- Output type → any input type (required or optional) ✓

**Record Merge (`r1 | r2`):**
- Same field name: types must match exactly
- Right record overrides left (feature)

**Type Matrix:**

If functions `a` and `b` are pipelined together,
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

### 4. Workflow Input Inference

- Workflow input type is inferred from the unmet inputs of all tasks
- An input for a task is met if there is an upstream task with an output by the
  same name

---

## Semantics

### Records
- Collection of key-value pairs
- Keys are field names (NAME)
- Values are strings or references

### Function Application
- `task record` applies task to record
- Partial application returns new function

### Pipeline
- `A |> B |> C` desugars to:
```
\x ->
    a = A x
    b = B (x | a)
    c = C (x | a | b)
    a | b | c
```

### Record Merge
- `r1 | r2` merges records, right overrides left
- Extra fields are allowed (ignored)

---

## Examples

### Minimal Workflow
```swl
align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

align |> sort |> call
```

### With Partial Application
```swl
align = import "align.sh"

align_hg38 = align {
    ref: "hg38.fa",
    ref_ann: "hg38.fa.ann",
}

align_hg38 { fastq1: "sample1.r1.fq", fastq2: "sample.r2.fq" }
```

### Explicit Arguments
```swl
\x ->
    a = align x
    s = sort { bam: a.bam, outbase: x.outbase }
    c = call { bam: s.bam, ref: x.ref, outbase: x.outbase }

    { bam: s.bam, bai: s.bai, bcf: c.bcf }
```

