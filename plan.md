# Implementation Plan

## Phase 1: Lexer - COMPLETE

### 1.1 Update Workflow Lexer - DONE
- Add `update` token type for `//` (record update)
- Change `pipe` token to `chain` for `|` (function chain)

### 1.2 Create Annotation Lexer - DONE (part of annotation.py)
- Parse bash scripts, extract comment lines
- Extract: task-doc, sections (in/out/run), params

## Phase 2: Parser - COMPLETE

### 2.1 Update Workflow Parser - DONE
- Parse `//` as record-update operator
- Parse `|` as function-chain operator
- Update precedence:
  1. Application (ws)
  2. Record update (`//`)
  3. Lambda arrow (`->`)
  4. Function chain (`|`)
  5. Let binding (`=`)

### 2.2 Create Annotation Parser - DONE
- Parse extracted comments into structured annotations

## Phase 3: AST - COMPLETE

### 3.1 Update Node Types - DONE
- Update for `//` 
- Chain for `|`
- Import for built-in import

### 3.2 Annotation AST Types - DONE
- Annotation: task-doc + sections
- Section: in/out/run + params
- Parameter: names, type, default, description

## Phase 4: Type System - COMPLETE

### 4.1 Type Checking - DONE
- Type compatibility in chains
- Type compatibility in unions
- Workflow input inference

### 4.2 Type Representation - DONE
- Simple types: `file`, `str`, `int`, `float`
- Optional: `file?`, `str?`, `int?`, `float?`
- Array types: `[file]`, `[str]`, `[int]`, `[float]`

---

## Remaining Work

### Phase 5: Interpreter (not started)
- Implement import (load task/workflow file, parse annotations, return function)
- Evaluate AST to workflow DAG
- Execute tasks in dependency order
- Handle Docker images from annotations
