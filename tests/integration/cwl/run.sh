#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
INT_DIR="$ROOT/tests/integration/cwl"
SWL_DIR="$ROOT/tests/integration/swl"
JOBS_DIR="$INT_DIR/jobs"
INPUTS_DIR="$ROOT/tests/integration/inputs"
OUTPUTS_DIR="$INT_DIR/outputs"
PASS=0
FAIL=0

cwltool --version

mkdir -p "$OUTPUTS_DIR"

export PYTHONPATH="$ROOT/python"
for swl in function pipe explicit panel map map_by; do
  python -m swl.compile "$SWL_DIR/${swl}.swl" -o "$INT_DIR/${swl}.json"
done

for swl in function pipe explicit panel map map_by; do
  python -m swl.transpile.cwl "$INT_DIR/${swl}.json" -o "$INT_DIR/${swl}.cwl"
done

run_test() {
  local name="$1"
  local cwl="$INT_DIR/${name}.cwl"
  local job="$JOBS_DIR/${name}.yml"
  local outdir="$OUTPUTS_DIR/$name"
  mkdir -p "$outdir"
  echo "=== $name ==="
  if cwltool --outdir "$outdir" "$cwl" "$job" 2>"$outdir/stderr.log"; then
    echo "PASS: $name"
    PASS=$((PASS + 1))
  else
    echo "FAIL: $name (see $outdir/stderr.log)"
    FAIL=$((FAIL + 1))
  fi
}

run_test function
run_test panel
run_test map
run_test map_by

echo "---"
echo "Passed: $PASS / $((PASS + FAIL))"
exit $FAIL
