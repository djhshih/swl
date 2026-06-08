#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
INT_DIR="$ROOT/tests/integration/cwl"
SWL_DIR="$ROOT/tests/integration/swl"
JOBS_DIR="$INT_DIR/jobs"
INPUTS_DIR="$INT_DIR/inputs"
OUTPUTS_DIR="$INT_DIR/outputs"
PASS=0
FAIL=0

cwltool --version

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR"
for f in sample1.r1.fq sample1.r2.fq sample2.r1.fq sample2.r2.fq \
         ref.fa ref.fa.amb ref.fa.ann ref.fa.bwt \
         ref.fa.fai ref.fa.pac ref.fa.sa empty.fq; do
  echo "$f" > "$INPUTS_DIR/$f"
done

export PYTHONPATH="$ROOT/python"
for swl in function pipe explicit panel map; do
  python -m swl.compile "$SWL_DIR/${swl}.swl" -o "$INT_DIR/${swl}.json"
done

for swl in function pipe explicit panel map; do
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

echo "---"
echo "Passed: $PASS / $((PASS + FAIL))"
exit $FAIL
