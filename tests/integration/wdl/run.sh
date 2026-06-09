#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
INT_DIR="$ROOT/tests/integration/wdl"
SWL_DIR="$ROOT/tests/integration/swl"
JOBS_DIR="$INT_DIR/jobs"
INPUTS_DIR="$ROOT/tests/integration/inputs"
OUTPUTS_DIR="$INT_DIR/outputs"
PASS=0
FAIL=0

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$JOBS_DIR"

cromwell --version 2>&1 || true
womtool --version 2>&1 || true

export PYTHONPATH="$ROOT/python"

# compile SWL to DAG JSON
for swl in function panel map map_by; do
  python -m swl.compile "$SWL_DIR/${swl}.swl" -o "$INT_DIR/${swl}.json"
done

# transpile DAG to WDL
for swl in function panel map map_by; do
  python -m swl.transpile.wdl "$INT_DIR/${swl}.json" -o "$INT_DIR/${swl}.wdl"
done

cromwell_run() {
  local name="$1"
  local wdl="$INT_DIR/${name}.wdl"
  local json="$JOBS_DIR/${name}.json"
  local outdir="$OUTPUTS_DIR/$name"
  mkdir -p "$outdir"
  echo "=== $name ==="

  # validate WDL with womtool
  if ! womtool validate "$wdl" 2>"$outdir/womtool_stderr.log"; then
    echo "FAIL: $name (womtool validation failed)"
    FAIL=$((FAIL + 1))
    return
  fi

  # run with Cromwell
  local cromwell_out="$outdir/cromwell_out"
  mkdir -p "$cromwell_out"
  if (
    cd "$cromwell_out" && \
    cromwell run "$wdl" -i "$json" > "$outdir/stdout.log" 2>"$outdir/stderr.log"
  ); then
    echo "PASS: $name"
    PASS=$((PASS + 1))
  else
    echo "FAIL: $name (see $outdir/stderr.log)"
    FAIL=$((FAIL + 1))
  fi
}

# generate JSON input files with absolute paths
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }

S1R1=$(abs "$INPUTS_DIR/sample1.r1.fq")
S1R2=$(abs "$INPUTS_DIR/sample1.r2.fq")
S2R1=$(abs "$INPUTS_DIR/sample2.r1.fq")
S2R2=$(abs "$INPUTS_DIR/sample2.r2.fq")
REF=$(abs "$INPUTS_DIR/ref.fa")
REF_AMB=$(abs "$INPUTS_DIR/ref.fa.amb")
REF_ANN=$(abs "$INPUTS_DIR/ref.fa.ann")
REF_BWT=$(abs "$INPUTS_DIR/ref.fa.bwt")
REF_FAI=$(abs "$INPUTS_DIR/ref.fa.fai")
REF_PAC=$(abs "$INPUTS_DIR/ref.fa.pac")
REF_SA=$(abs "$INPUTS_DIR/ref.fa.sa")
EMPTY=$(abs "$INPUTS_DIR/empty.fq")

# WDL expects plain file paths (not CWL-style {"class":"File","path":"..."})
# --- function.json ---
cat > "$JOBS_DIR/function.json" << EOF
{
  "main.fastq1": "$S1R1",
  "main.fastq2": "$S1R2",
  "main.outbase": "test",
  "main.ref": "$REF",
  "main.ref_amb": "$REF_AMB",
  "main.ref_ann": "$REF_ANN",
  "main.ref_bwt": "$REF_BWT",
  "main.ref_fai": "$REF_FAI",
  "main.ref_pac": "$REF_PAC",
  "main.ref_sa": "$REF_SA"
}
EOF

# panel and map use sub-workflows which Cromwell 86 does not support;
# map_by requires collect_by_key (WDL 1.1). Only function is tested end-to-end.

cromwell_run function

echo "---"
echo "Passed: $PASS / $((PASS + FAIL))"
exit $FAIL
