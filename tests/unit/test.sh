#!/bin/bash
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
echo "running unit tests"
PYTHONPATH="$ROOT/python:$ROOT/tests/unit" python -m unittest \
    tests.unit.swl.syntax.wf.test_lexer \
    tests.unit.swl.syntax.wf.test_parser \
    tests.unit.swl.syntax.task.test_parser \
    tests.unit.swl.syntax.task.test_interpolation \
    tests.unit.swl.syntax.task.test_bash \
    tests.unit.swl.semantic.task.test_type \
    tests.unit.swl.semantic.wf.test_check \
    tests.unit.swl.ir.test_lower \
    tests.unit.swl.ir.test_force \
    tests.unit.swl.ir.test_force_codec \
    tests.unit.swl.transpile.cwl.test_emit \
    tests.unit.swl.transpile.wdl.test_emit \
    tests.unit.swl.transpile.nf.test_emit
printf "unit tests passed\n\n"
