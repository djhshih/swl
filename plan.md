# Implementation Plan

## Phase 1: Lexer

### 1.1 Update Workflow Lexer
- Add `update` token type for `//` (record update)
- Change `pipe` token to `chain` for `|` (function chain)

### 1.2 Create Annotation Lexer
- Parse bash scripts, extract comment lines
- Strip `#` prefix per spec assumption
- Extract: task-doc, sections (in/out/run), params

## Phase 2: Parser

### 2.1 Update Workflow Parser
- Parse `//` as record-update operator
- Parse `|` as function-chain operator
- Update precedence:
  1. Application (ws)
  2. Record update (`//`)
  3. Lambda arrow (`->`)
  4. Function chain (`|`)
  5. Let binding (`=`)

### 2.2 Create Annotation Parser
- Parse extracted comments into structured annotations

## Phase 3: AST

### 3.1 Update Node Types
- RecordUpdate for `//` (was Update)
- FunctionChain for `|` (was Pipe)
- Import (for built-in import function)

### 3.2 Annotation AST Types
- Annotation: task-doc + sections
- Section: in/out/run + params
- Parameter: names, type, default, description

## Phase 4: Type System

### 4.1 Type Checking
- Implement compile-time checks from spec:
  1. Import verification (file exists, no cycles)
  2. Type compatibility in chains/unions
  3. DAG circularity check
  4. Workflow input inference

### 4.2 Type Representation
- Simple types: `file`, `str`, `int`, `float`
- Optional: `file?`, `str?`, `int?`, `float?`
- Array types: `[file]`, `[str]`, `[int]`, `[float]`

## Phase 5: Interpreter

### 5.1 Built-in Functions
- `import` - loads task/workflow file, returns function

### 5.2 Execution
- Evaluate AST to workflow DAG
- Execute tasks in dependency order
- Handle Docker images from annotations

---

## Files to Modify

| File | Action | Description |
|------|--------|-------------|
| lexer.py | Modify | Add `//` token, update `|` token |
| parser.py | Modify | Update operators, precedence |
| node.py | Modify | Rename Update→RecordUpdate, Pipe→Chain |
| eval.py | Modify | Add import, type checking |

## Files to Create

| File | Description |
|------|-------------|
| annotation.py | Parse task annotations |
| type.py | Type representations and checking |
