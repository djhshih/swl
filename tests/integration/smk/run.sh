#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
INT_DIR="$ROOT/tests/integration/smk"
DAG_DIR="$ROOT/tests/integration/dag"
CONFIG_DIR="$INT_DIR/config"
OUTPUTS_DIR="$INT_DIR/outputs"
PASS=0
FAIL=0

if command -v snakemake &>/dev/null; then
    SNAKEMAKE="snakemake"
    snakemake --version 2>&1 | head -1
fi

mkdir -p "$OUTPUTS_DIR"
export PYTHONPATH="$ROOT/python"

for swl in function panel map map_by; do
    python -m swl.transpile.smk "$DAG_DIR/${swl}.json" -o "$INT_DIR/${swl}.smk"
done

run_test() {
    local name="$1"
    local smk="$INT_DIR/${name}.smk"
    local config="$CONFIG_DIR/${name}.yaml"
    local outdir="$OUTPUTS_DIR/$name"
    mkdir -p "$outdir"
    local log="$outdir/snakemake.log"

    echo "=== $name ==="

    if $SNAKEMAKE -nq -s "$smk" --configfile "$config" --directory "$INT_DIR" > "$log" 2>&1; then
        echo "LINT PASS: $name (lint)"
        PASS=$((PASS + 1))
    elif grep -q "SyntaxError" "$log"; then
        echo "LINT FAIL: $name (syntax error, see $log)"
        FAIL=$((FAIL + 1))
    else
        echo "LINT FAIL: $name (see $log)"
        FAIL=$((FAIL + 1))
    fi

    if $SNAKEMAKE -nq -s "$smk" --configfile "$config" --directory "$INT_DIR" > "$log" 2>&1; then
        echo "PASS: $name"
        PASS=$((PASS + 1))
    elif grep -q "KeyError" "$log"; then
        echo "FAIL: $name (missing config key, see $log)"
        FAIL=$((FAIL + 1))
    else
        echo "FAIL: $name (see $log)"
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
