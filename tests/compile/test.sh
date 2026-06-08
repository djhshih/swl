#!/bin/bash
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT/tests/compile/swl"

TASK_COUNT=0
WF_COUNT=0
COMPARE_COUNT=0
CWL_COMPARE_COUNT=0
WDL_COMPARE_COUNT=0
NF_COMPARE_COUNT=0
SMK_COMPARE_COUNT=0

EXPECT_COMPILE_FAIL=(
    bad_explicit.swl
    bad_pipe.swl
)

expects_compile_fail() {
    local target="$1"
    for item in "${EXPECT_COMPILE_FAIL[@]}"; do
        if [[ "$item" == "$target" ]]; then return 0; fi
    done
    return 1
}

evaluate_task() {
    TASK_COUNT=$((TASK_COUNT + 1))
    set +e
    echo "file: $1"
    echo "syntax:"
    PYTHONPATH="$ROOT/python" python -m swl.eval.syntax_task "$1"
    printf "return: $?\n"
    echo "semantic:"
    PYTHONPATH="$ROOT/python" python -m swl.eval.semantic_task "$1"
    printf "return: $?\n\n"
    set -e
}

evaluate_wf() {
    WF_COUNT=$((WF_COUNT + 1))
    set +e
    local wf="$1"
    local out="../dag/$(basename "${wf%.swl}").json"
    local expect_fail=0
    if expects_compile_fail "$wf"; then expect_fail=1; fi

    echo "file: $wf"
    if (( expect_fail )); then echo "expectation: compile failure"
    else echo "expectation: compile success"; fi

    echo "syntax:"
    PYTHONPATH="$ROOT/python" python -m swl.eval.syntax_wf "$wf"
    printf "return: $?\n"
    echo "semantic:"
    PYTHONPATH="$ROOT/python" python -m swl.eval.semantic_wf "$wf"
    printf "return: $?\n"
    echo "ir:"
    PYTHONPATH="$ROOT/python" python -m swl.eval.ir "$wf"
    local ir_status=$?
    printf "return: %d" "$ir_status"
    if (( expect_fail )) && (( ir_status != 0 )); then printf " (expected failure)"; fi
    printf "\n"
    echo "dag:"
    PYTHONPATH="$ROOT/python" python -m swl.eval.dag "$wf"
    local dag_status=$?
    printf "return: %d" "$dag_status"
    if (( expect_fail )) && (( dag_status != 0 )); then printf " (expected failure)"; fi
    printf "\n"
    echo "compile:"
    rm -f "$out"
    PYTHONPATH="$ROOT/python" python -m swl.compile "$wf" -o "$out"
    local compile_status=$?
    printf "return: %d" "$compile_status"
    if (( expect_fail )); then
        if (( compile_status != 0 )); then printf " (expected failure)"
        else printf " (UNEXPECTED SUCCESS)"; set -e; return 1; fi
        if [[ -e "$out" ]]; then echo; echo "unexpected artifact: $out"; set -e; return 1; fi
    else
        if (( compile_status == 0 )); then printf " (expected success)"
        else printf " (UNEXPECTED FAILURE)"; set -e; return 1; fi
        if [[ ! -e "$out" ]]; then echo; echo "missing artifact: $out"; set -e; return 1; fi
    fi
    printf "\n\n"
    set -e
}

compare_dags() {
    local left="$1" right="$2"
    COMPARE_COUNT=$((COMPARE_COUNT + 1))
    echo "compare: $left == $right"
    PYTHONPATH="$ROOT/python" python - <<'PY' "$left" "$right"
import json, sys
left = json.load(open(sys.argv[1]))
right = json.load(open(sys.argv[2]))
if left != right: raise SystemExit(1)
PY
    echo "result: ok"; echo
}

compare_files() {
    local left="$1" right="$2"
    COMPARE_COUNT=$((COMPARE_COUNT + 1))
    echo "compare: $left == $right"
    diff -u "$left" "$right" >/dev/null
    echo "result: ok"; echo
}

print_summary() {
    echo "summary:"
    echo "  task files checked: $TASK_COUNT"
    echo "  workflow files checked: $WF_COUNT"
    echo "  dag equality checks: $COMPARE_COUNT"
    echo "  cwl golden checks: $CWL_COMPARE_COUNT"
    echo "  wdl golden checks: $WDL_COMPARE_COUNT"
    echo "  nf golden checks: $NF_COMPARE_COUNT"
    echo "  smk golden checks: $SMK_COMPARE_COUNT"
}

if (( $# > 0 )); then
    if [[ $1 =~ \.swl ]]; then evaluate_wf "$1"
    else evaluate_task "$1"; fi
else
    mkdir -p ../dag && rm -f ../dag/*
    for task in *.sh; do evaluate_task "$task"; done
    for swl in *.swl; do evaluate_wf "$swl"; done

    compare_dags ../dag/pipe.json ../dag/function.json
    compare_dags ../dag/explicit.json ../dag/function.json

    mkdir -p ../cwl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.cwl ../dag/function.json -o ../cwl/function.cwl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.cwl ../dag/pipe.json    -o ../cwl/pipe.cwl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.cwl ../dag/explicit.json -o ../cwl/explicit.cwl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.cwl ../dag/panel.json   -o ../cwl/panel.cwl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.cwl ../dag/map.json     -o ../cwl/map.cwl

    compare_files ../cwl/pipe.cwl ../cwl/function.cwl
    compare_files ../cwl/explicit.cwl ../cwl/function.cwl
    CWL_COMPARE_COUNT=$((CWL_COMPARE_COUNT + 2))

    mkdir -p ../wdl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.wdl ../dag/function.json   -o ../wdl/function.wdl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.wdl ../dag/pipe.json       -o ../wdl/pipe.wdl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.wdl ../dag/explicit.json   -o ../wdl/explicit.wdl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.wdl ../dag/panel.json      -o ../wdl/panel.wdl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.wdl ../dag/map.json        -o ../wdl/map.wdl
    PYTHONPATH="$ROOT/python" python -m swl.transpile.wdl ../dag/map_by.json     -o ../wdl/map_by.wdl

    compare_files ../wdl/pipe.wdl ../wdl/function.wdl
    compare_files ../wdl/explicit.wdl ../wdl/function.wdl
    WDL_COMPARE_COUNT=$((WDL_COMPARE_COUNT + 3))

    mkdir -p ../nf
    PYTHONPATH="$ROOT/python" python -m swl.transpile.nf ../dag/function.json    -o ../nf/function.nf
    PYTHONPATH="$ROOT/python" python -m swl.transpile.nf ../dag/pipe.json        -o ../nf/pipe.nf
    PYTHONPATH="$ROOT/python" python -m swl.transpile.nf ../dag/explicit.json    -o ../nf/explicit.nf
    PYTHONPATH="$ROOT/python" python -m swl.transpile.nf ../dag/panel.json       -o ../nf/panel.nf
    PYTHONPATH="$ROOT/python" python -m swl.transpile.nf ../dag/map.json         -o ../nf/map.nf
    PYTHONPATH="$ROOT/python" python -m swl.transpile.nf ../dag/map_by.json      -o ../nf/map_by.nf

    compare_files ../nf/pipe.nf ../nf/function.nf
    compare_files ../nf/explicit.nf ../nf/function.nf
    NF_COMPARE_COUNT=$((NF_COMPARE_COUNT + 2))

    mkdir -p ../smk
    PYTHONPATH="$ROOT/python" python -m swl.transpile.smk ../dag/function.json   -o ../smk/function.smk
    PYTHONPATH="$ROOT/python" python -m swl.transpile.smk ../dag/pipe.json       -o ../smk/pipe.smk
    PYTHONPATH="$ROOT/python" python -m swl.transpile.smk ../dag/explicit.json   -o ../smk/explicit.smk
    PYTHONPATH="$ROOT/python" python -m swl.transpile.smk ../dag/panel.json      -o ../smk/panel.smk
    PYTHONPATH="$ROOT/python" python -m swl.transpile.smk ../dag/map.json        -o ../smk/map.smk
    PYTHONPATH="$ROOT/python" python -m swl.transpile.smk ../dag/map_by.json     -o ../smk/map_by.smk

    compare_files ../smk/pipe.smk ../smk/function.smk
    compare_files ../smk/explicit.smk ../smk/function.smk
    SMK_COMPARE_COUNT=$((SMK_COMPARE_COUNT + 2))

    print_summary
fi
