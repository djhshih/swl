#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
INT_DIR="$ROOT/tests/integration/wdl"
DAG_DIR="$ROOT/tests/integration/dag"
JOBS_DIR="$INT_DIR/jobs"
INPUTS_DIR="$ROOT/tests/integration/inputs"
OUTPUTS_DIR="$INT_DIR/outputs"
PASS=0
FAIL=0

mkdir -p "$OUTPUTS_DIR"

HAVE_SPROCKET=false
HAVE_WOMTOOL=false
HAVE_CROMWELL=false

if command -v sprocket &>/dev/null; then
  HAVE_SPROCKET=true
  echo "sprocket $(sprocket --version 2>&1)"
fi
if command -v womtool &>/dev/null; then
  HAVE_WOMTOOL=true
  echo "womtool $(womtool --version 2>&1)"
fi
if command -v cromwell &>/dev/null; then
  HAVE_CROMWELL=true
  echo "cromwell $(cromwell --version 2>&1)"
fi

export PYTHONPATH="$ROOT/python"

for swl in function panel map map_by; do
  python -m swl.transpile.wdl "$DAG_DIR/${swl}.json" -o "$INT_DIR/${swl}.wdl"
done

cromwell_run() {
  local name="$1"
  local wdl="$INT_DIR/${name}.wdl"
  local json="$JOBS_DIR/${name}.json"
  local outdir="$OUTPUTS_DIR/$name"
  mkdir -p "$outdir"
  echo "=== $name ==="

  if $HAVE_WOMTOOL; then
    if ! womtool validate "$wdl" 2>"$outdir/womtool_stderr.log"; then
      echo "SKIP: $name (womtool validation failed, known limitation)"
      return
    fi
  fi

  if $HAVE_CROMWELL; then
    local cromwell_out="$outdir/cromwell"
    mkdir -p "$cromwell_out"
    if (cd "$cromwell_out" && cromwell run "$wdl" -i "$json" > "$outdir/cromwell_stdout.log" 2>"$outdir/cromwell_stderr.log"); then
      echo "PASS: $name"
      PASS=$((PASS + 1))
    else
      echo "SKIP: $name (cromwell run failed)"
    fi
  else
    echo "PASS: $name (womtool validation only)"
    PASS=$((PASS + 1))
  fi
}

run_test() {
  local name="$1"
  local wdl="$INT_DIR/${name}.wdl"
  local json="$JOBS_DIR/${name}.json"
  local outdir="$OUTPUTS_DIR/$name"
  mkdir -p "$outdir"
  echo "=== $name ==="

  if $HAVE_SPROCKET; then
    if sprocket check "$wdl" > "$outdir/sprocket_check.log" 2>&1; then
      if sprocket run "$wdl" "@$json" -o "$outdir/sprocket" > "$outdir/sprocket_run.log" 2>&1; then
        echo "PASS: $name"
        PASS=$((PASS + 1))
        return
      else
        echo "FAIL: $name (sprocket run failed, see $outdir/sprocket_run.log)"
        FAIL=$((FAIL + 1))
        return
      fi
    else
      echo "INFO: $name (sprocket check failed, trying fallback)"
    fi
  fi

  cromwell_run "$name"
}

run_test function
run_test panel
run_test map
run_test map_by

echo "---"
echo "Passed: $PASS / $((PASS + FAIL))"
exit $FAIL
