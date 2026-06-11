# Test Improvement Plan

## What compile & integration tests cover

`tests/compile/test.sh` runs every `.swl`/`.sh` file through `eval.syntax_task`, `eval.semantic_task`, `eval.syntax_wf`, `eval.semantic_wf`, `eval.ir`, `eval.dag`, and `swl.compile`. It also transpiles every compiled DAG through all 4 backends (CWL/WDL/NF/SMK) and verifies pipe==function and explicit==function equivalence.

**Already covered — do NOT add redundant unit tests:**
- `swl/eval/*.py` — run on every SWL file
- `swl/compile.py` — CLI exercised directly
- All 4 transpile backends — end-to-end equivalence verified for 6 workflows
- All transpile CLI/`__main__.py` — invoked by compile test
- `swl.api.compile_workflow` — called via CLI
- `swl/dag/merge.py` — all functions private (`_`), exercised by force tests
- `swl/dag/tooldefs.py` — all functions private, exercised by force tests
- `swl/semantic/wf/bashvars.py` — private function, exercised via tooldefs
- `swl/semantic/wf/infer.py` — exercised by `test_check.py` (which calls `apply_function`, `application_result`, etc.)

---

## 1. Make Restrictive Tests Less Restrictive

These tests currently impede refactoring by checking deep implementation structure rather than behavioral contracts. Each change preserves the test's original intent while removing coupling to internal details.

### 1.1 `tests/unit/swl/ir/test_force.py`

**`test_force_saturated_workflow_produces_task_dag`** (line 266)
- Replace hardcoded step IDs `['align', 'sort']` with step count + output key check:

```python
self.assertEqual(len(data['steps']), 2)
self.assertEqual(sorted(data['outputs'].keys()), ['bai', 'bam'])
```

**`test_chain_root_is_instantiated_during_force`** (line 302, 16 assertions)
- Cut to 3-4: step count == 3, output count == 3, one dependency check:

```python
self.assertEqual(len(data['steps']), 3)
self.assertEqual(len(data['outputs']), 3)
self.assertEqual(data['steps'][2]['deps'], [data['steps'][1]['id']])
```

**`test_serialized_dag_is_self_contained`** (line 320)
- Replace exact path assertion with `assertTrue(data['steps'][0]['path'].endswith('align.sh'))`.

**`test_function_and_chain_compile_to_same_shape`** (line 331)
- Replace tuple comparison with `assertEqual(len(chain['steps']), len(function['steps']))` + compare just the structural summary (keys sets, dep counts).

**`test_force_partial_map_root_materializes_batch_workflow`** (line 285)
- Check `'xs'` ∈ inputs, step count == 1, map source type == `'input'`. Drop exact key-list assertions.

**All dedup tests** (lines 342–382: `test_reused_variable_...`, `test_partial_application_reuse_...`, `test_workflow_partial_application_reuse_...`, `test_nested_workflow_value_reuse_...`)
- Replace `[task['id'] for task in data['steps']] == ['align']` with `self.assertEqual(len(data['steps']), 1)`.
- For `test_workflow_partial_application_reuse` (which checks `['align', 'align_2']`): check step count == 2, not exact names.

**`test_mapped_table_source_uses_explicit_logical_table_metadata`** (line 444)
- Check source type is `'table'`, verify column count, spot-check one column — drop full dict equality.

### 1.2 `tests/unit/swl/ir/test_lower.py`

**`test_lower_imports_to_functions`** (line 93, 13 assertions into IR node internals)
- Cut to: verify Lambda, verify Block with 2 bindings.

```python
self.assertIsInstance(tree, ir.Lambda)
self.assertIsInstance(tree.body, ir.Block)
self.assertEqual(len(tree.body.bindings), 2)
```

**`test_lower_bindings_produce_variables_and_refs`** (line 113)
- Merge into `test_lower_imports_to_functions` (same input). Tests overlap.

### 1.3 `tests/unit/swl/transpile/cwl/test_emit.py`

**`test_transpile_function_workflow`** (line 157, 22 assertions)
- Check `cwlVersion`, count tools (== 3), count workflow outputs (== 3), verify one tool has `baseCommand` and `coresMin`:

```python
self.assertEqual(cwl['cwlVersion'], 'v1.2')
tools = [item for item in cwl['$graph'] if item['class'] == 'CommandLineTool']
self.assertEqual(len(tools), 3)
outputs = [item for item in cwl['$graph'][-1]['outputs']]
self.assertEqual(len(outputs), 3)
```

**`test_batch_mapped_task_emits_scatter_and_tab_column_input_type`** (line 249) + **`test_root_partial_map_transpiles_as_scattered_subworkflow`** (line 345)
- Replace 10-item scatter port lists with `self.assertIn('scatter', step)` + `self.assertGreater(len(step.get('scatter', [])), 0)`.

### 1.4 `tests/unit/swl/transpile/nf/test_emit.py`

**`test_transpile_function_workflow`** (line 152)
- Use `self.assertIn('process', nf)` instead of exact process-name assertions.

**`test_root_partial_map_inlines_pipeline`** (line 276)
- Replace `'SORT(ALIGN.out.bam, outbase)'` with `self.assertIn('ALIGN', nf)` + `self.assertIn('SORT', nf)`.

### 1.5 `tests/unit/swl/transpile/wdl/test_emit.py`

**`test_transpile_function_workflow`** (line 161, 27 assertions)
- Keep 8–10: `version 1.0`, 3 `task` sections, `workflow main`, one input declaration, one output declaration.

---

## 2. Improve Trivial Tests

### 2.1 `tests/unit/swl/syntax/task/test_interpolation.py`

Currently 5 tests. The `Parser.parse_word` method is a public API with contract-based output (returns `Word`). These additions close edge-case gaps:

```python
def test_literal_before_var(self):
    result = Parser().parse_word('prefix${var}')
    self.assertEqual(result, Word([Literal('prefix'), Var('var')]))

def test_multiple_vars(self):
    result = Parser().parse_word('${a}${b}')
    self.assertEqual(result, Word([Var('a'), Var('b')]))

def test_escaped_dollar_is_literal(self):
    result = Parser().parse_word('\\${notavar}')
    self.assertEqual(result, Word([Literal('${notavar}')]))

def test_empty_brace_fails(self):
    with self.assertRaises(ValueError):
        Parser().parse_word('${}')

def test_complex_expression(self):
    result = Parser().parse_word('${a + b * c}')
    self.assertEqual(result, Word([Expr('a + b * c')]))

def test_empty_string(self):
    result = Parser().parse_word('')
    self.assertEqual(result, Word([]))
```

### 2.2 `tests/unit/swl/syntax/wf/test_lexer.py`

Currently 7 tests, all happy-path token sequences. `Lexer` is a public class:

```python
def test_empty_string(self):
    lexer = Lexer('')
    self.assertEqual([x for x in lexer], [Token(TokenType.eof)])

def test_unterminated_string_fails(self):
    lexer = Lexer('x = "unterminated')
    with self.assertRaises(ValueError):
        [x for x in lexer]

def test_just_comment_line(self):
    lexer = Lexer('# just a comment')
    self.assertEqual([x for x in lexer], [Token(TokenType.eof)])

def test_multiple_blank_lines(self):
    lexer = Lexer('\n\n\n')
    self.assertEqual([x for x in lexer], [Token(TokenType.eof)])
```

### 2.3 `tests/unit/swl/syntax/task/test_bash.py`

Currently 2 tests. `bash.parse` is the public API:

```python
def test_empty_script(self):
    script = bash.parse('')
    self.assertEqual(len(script.statements), 0)

def test_script_with_only_comment(self):
    script = bash.parse('# just a comment\n')
    self.assertEqual(len(script.statements), 0)
```

### 2.4 `tests/unit/swl/syntax/wf/test_parser.py`

10 of 14 tests use `assertIsNotNone` — they verify parsing "doesn't crash" but never check the result. Replace with shallow top-level node type checks:

```python
from swl.syntax.wf.node import NodeType as NT

result = Parser().parse(src)
self.assertEqual(result.type, NT.fun)   # for lambda tests
# or
self.assertEqual(result.type, NT.id)    # for simple expression tests
# or
self.assertEqual(result.type, NT.chain) # for chain tests
```

Use the appropriate `NodeType` per test case (see `swl/syntax/wf/node.py` for enum values). These verify what shape the parser returns without descending into subtrees.

### 2.5 `tests/unit/swl/semantic/wf/test_validate.py`

Currently 4 tests. `validate_workflow_inputs` is a public function:

```python
def test_validate_missing_input_reports_error(self):
    root = self._fixture_dir()
    self._write(root, 'align.sh', _ALIGN)
    path = self._write(root, 'simple.swl', _SIMPLE)
    result = Checker().load(path)
    with self.assertRaises(WorkflowInputValidationError) as ctx:
        validate_workflow_inputs(result, {'outbase': 'x'})
    self.assertIn('fastq1', str(ctx.exception))

def test_validate_batch_single_element_arrays(self):
    root = self._fixture_dir()
    self._write(root, 'align.sh', _ALIGN)
    self._write(root, 'merge.sh', _MERGE)
    path = self._write(root, 'batch.swl', _BATCH)
    result = Checker().load(path)
    value = {'fastq1': ['a.fq'], 'fastq2': ['b.fq'], 'ref': ['hg38.fa'],
             'ref_fai': ['hg38.fa.fai'], 'outbase': ['a']}
    self.assertEqual(validate_workflow_inputs(result, value), value)
```

---

## 3. Fill Coverage Gaps

Only add tests for code genuinely not exercised by compile/integration tests.

### 3.1 `tests/unit/swl/test_api.py` (NEW)

Three API functions are never called by compile tests (they call the underlying modules directly):

```python
import unittest, json, os, tempfile
from swl.api import force_workflow, load_workflow, transpile_dag
from swl.dag.node import DAG

class TestAPI(unittest.TestCase):
    def _files(self):
        return {
            '/v/align.sh': '# @ Align\n# in\n#   x file\n# out\n#   y file = out.txt\necho hi\n',
            '/v/wf.swl': 'a = import "align.sh"\na\n',
        }

    def test_force_workflow_returns_dag(self):
        dag = force_workflow('/v/wf.swl', files=self._files())
        self.assertIsInstance(dag, DAG)
        self.assertEqual(len(dag.steps), 1)

    def test_force_workflow_errors_on_invalid(self):
        with self.assertRaises(Exception):
            force_workflow('/v/nonexistent.swl', files=self._files())

    def test_load_workflow_returns_check_result(self):
        result = load_workflow('/v/wf.swl', files=self._files())
        self.assertIsNotNone(result.signature)
        self.assertIn('y', result.signature.outputs)

    def test_transpile_dag_cwl(self):
        dag = force_workflow('/v/wf.swl', files=self._files())
        td = tempfile.TemporaryDirectory()
        dag.write(os.path.join(td.name, 'plan.json'))
        result = transpile_dag(os.path.join(td.name, 'plan.json'), 'cwl')
        parsed = json.loads(result)
        self.assertEqual(parsed['cwlVersion'], 'v1.2')

    def test_transpile_dag_nf(self):
        dag = force_workflow('/v/wf.swl', files=self._files())
        td = tempfile.TemporaryDirectory()
        dag.write(os.path.join(td.name, 'plan.json'))
        result = transpile_dag(os.path.join(td.name, 'plan.json'), 'nf')
        self.assertIn('process', result)

    def test_transpile_dag_invalid_target(self):
        with self.assertRaises(ValueError):
            transpile_dag('/fake.json', 'invalid')
```

### 3.2 Expand test runner

Add new module to `tests/unit/test.sh`:

```bash
tests.unit.swl.test_api \
```

---

## 4. Cleanup

- `tests/unit/swl/semantic/wf/test_check.py` line 663: Remove commented-out dead test:
  ```python
  # def test_table_update_reports_explicit_not_implemented(self):
  ```

---

## Summary

| Category | Change | Refactoring impact |
|----------|--------|--------------------|
| **Harden** 15 brittle tests | Replace deep assertion chains with count/type/contract checks | Tests break only when external contract changes, not when internals are reorganized |
| **Deepen** 5 trivial test files | Add edge-case coverage (empty input, error paths, boundary conditions) | Tests exercise public APIs with inputs the compile tests never use |
| **New** `test_api.py` | Test 3 truly uncovered public API functions | Catches regressions in the public interface |
| **Remove** dead code | Delete commented-out test | N/A |

**Deliberately excluded** (covered by compile/integration or impedes refactoring):
- `merge.py`, `tooldefs.py`, `bashvars.py` — all private functions, exercised by force/compile tests
- `infer.py` — exercised by `test_check.py`
- `loader.py` — exercised by compile (Checker.load uses Loader)
- All `eval/*` modules — run on every SWL file in compile tests
- `common.py` — exercised by all 4 transpile backends across 6 workflows each
- All transpile backends — end-to-end equivalence verified by compile tests
- All CLI/`__main__.py` — invoked by compile test runner
