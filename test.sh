#!/bin/bash
set -e

PASS=0
FAIL=0

run_suite() {
    local name="$1" path="$2"
    echo "===== $name ====="
    if bash "$path"; then
        echo "PASS: $name"
        PASS=$((PASS + 1))
    else
        echo "FAIL: $name"
        FAIL=$((FAIL + 1))
    fi
    echo
}

run_suite "unit"       "tests/unit/test.sh"
run_suite "compile"    "tests/compile/test.sh"
run_suite "integration-cwl" "tests/integration/cwl/run.sh"
run_suite "integration-wdl" "tests/integration/wdl/run.sh"

echo "===== RESULTS ====="
echo "Passed: $PASS / $((PASS + FAIL))"
