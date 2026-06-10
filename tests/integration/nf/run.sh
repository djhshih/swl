#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
INT_DIR="$ROOT/tests/integration/nf"
DAG_DIR="$ROOT/tests/integration/dag"
PARAMS_DIR="$INT_DIR/params"
OUTPUTS_DIR="$INT_DIR/outputs"
PASS=0
FAIL=0

NEXTFLOW=""
if command -v nextflow &>/dev/null; then
    NEXTFLOW="nextflow"
    nextflow -version 2>&1 | head -2
fi

# Create dummy input files for actual nextflow execution
mkdir -p /tmp/dummy
for sample in sample1 sample2; do
    touch "/tmp/dummy/${sample}.r1.fq" "/tmp/dummy/${sample}.r2.fq"
done
touch /tmp/dummy/ref.fa{,.amb,.ann,.bwt,.fai,.pac,.sa}

mkdir -p "$OUTPUTS_DIR"
export PYTHONPATH="$ROOT/python"

for swl in function pipe explicit panel map map_by; do
    python -m swl.transpile.nf "$DAG_DIR/${swl}.json" -o "$INT_DIR/${swl}.nf"
done

run_test() {
    local name="$1"
    local nf="$INT_DIR/${name}.nf"
    local params="$PARAMS_DIR/${name}.json"
    local outdir="$OUTPUTS_DIR/$name"
    mkdir -p "$outdir"
    echo "=== $name ==="

    local lint_log="$outdir/lint.log"
    if nextflow lint "$nf" > "$lint_log" 2>&1; then
        echo "LINT PASS: $name"
    elif grep -q "error" "$lint_log"; then
        echo "LINT FAIL: $name (see $lint_log)"
        FAIL=$((FAIL + 1))
        return
    else
        echo "LINT PASS: $name"
    fi

    local run_log="$outdir/run.log"
    if $NEXTFLOW run "$nf" -params-file "$params" > "$run_log" 2>&1; then
        echo "PASS: $name"
        PASS=$((PASS + 1))
    else
        echo "FAIL: $name (run failed, see $run_log)"
        FAIL=$((FAIL + 1))
    fi
}

run_test function
run_test pipe
run_test explicit
run_test map
run_test panel
run_test map_by

echo "---"
echo "Passed: $PASS / $((PASS + FAIL))"
exit $FAIL
