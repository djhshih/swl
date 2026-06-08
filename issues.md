# Spec/Implementation Discrepancies

## 1. Missing `bool` type from SWL type system

**File:** `doc/spec.md:214-218` vs `doc/dag.md:220` vs `python/swl/dag/finalize.py:126`

The SWL spec defines only `file`, `str`, `int`, `float` as valid types. The DAG spec says `OutputSpec.type` accepts `"bool"`, `"bool?"`. The compiler emits `'bool'` from `_literal_type()` (finalize.py:126) when an output is a Python `bool`. The type system has no `bool` type — this is an undeclared extension.

