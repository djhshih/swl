#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SWL_DIR="$ROOT/tests/integration/swl"
DAG_DIR="$ROOT/tests/integration/dag"
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

export PYTHONPATH="$ROOT/python"
mkdir -p "$DAG_DIR"

for swl in function pipe explicit panel map map_by; do
  python -m swl.compile "$SWL_DIR/${swl}.swl" -o "$DAG_DIR/${swl}.json"
done

run_suite "integration-cwl" "$ROOT/tests/integration/cwl/run.sh"
run_suite "integration-wdl" "$ROOT/tests/integration/wdl/run.sh"

echo "===== RESULTS ====="
echo "Passed: $PASS / $((PASS + FAIL))"
exit $FAIL
