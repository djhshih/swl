# Implementation Plan

## Phase 1: Lexer - COMPLETE

### 1.2 Create Annotation Lexer - DONE (part of annotation.py)
- Parse bash scripts, extract comment lines
- Extract: task-doc, sections (in/out/run), params

### 2.2 Create Annotation Parser - DONE
- Parse extracted comments into structured annotations

## Phase 4: Type System

### 4.1 Type Checking
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
