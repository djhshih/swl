#!/bin/bash

set -e

TASK_COUNT=0
WF_COUNT=0
COMPARE_COUNT=0

EXPECT_COMPILE_FAIL=(
	tests/bad_explicit.swl
	tests/bad_pipe.swl
)

expects_compile_fail() {
	local target="$1"
	for item in "${EXPECT_COMPILE_FAIL[@]}"; do
		if [[ "$item" == "$target" ]]; then
			return 0
		fi
	done
	return 1
}

run_unit_tests() {
	echo "running unit tests"
	PYTHONPATH=python python -m unittest \
		swl.syntax.lexer \
		swl.syntax.parser \
		swl.syntax.task.test_parser \
		swl.syntax.task.test_interpolation \
		swl.syntax.task.test_bash \
		swl.semantic.task.test_type \
		swl.semantic.wf.test_check \
		swl.ir.test_lower \
		swl.ir.test_force \
		swl.ir.test_force_codec \
		swl.transpile.cwl.test_emit
	printf "unit tests passed\n\n"
}

evaluate_task() {
	TASK_COUNT=$((TASK_COUNT + 1))
	set +e
	echo "file: $1"
	echo "syntax:"
	PYTHONPATH=python python -m swl.eval_task $1
	printf "return: $?\n"
	echo "semantic:"
	PYTHONPATH=python python -m swl.eval_task_semantic $1
	printf "return: $?\n\n"
	set -e
}

evaluate_wf() {
	WF_COUNT=$((WF_COUNT + 1))
	set +e
	local wf="$1"
	local out="tests/dag/$(basename ${wf%.swl}).json"
	local expect_fail=0
	if expects_compile_fail "$wf"; then
		expect_fail=1
	fi

	echo "file: $wf"
	if (( expect_fail )); then
		echo "expectation: compile failure"
	else
		echo "expectation: compile success"
	fi

	echo "syntax:"
	PYTHONPATH=python python -m swl.eval "$wf"
	printf "return: $?\n"

	echo "semantic:"
	PYTHONPATH=python python -m swl.eval_wf_semantic "$wf"
	printf "return: $?\n"

	echo "ir:"
	PYTHONPATH=python python -m swl.eval_ir "$wf"
	local ir_status=$?
	printf "return: %d" "$ir_status"
	if (( expect_fail )) && (( ir_status != 0 )); then
		printf " (expected failure)"
	fi
	printf "\n"

	echo "force:"
	PYTHONPATH=python python -m swl.eval_force "$wf"
	local force_status=$?
	printf "return: %d" "$force_status"
	if (( expect_fail )) && (( force_status != 0 )); then
		printf " (expected failure)"
	fi
	printf "\n"

	echo "compile:"
	rm -f "$out"
	PYTHONPATH=python python -m swl.compile "$wf" -o "$out"
	local compile_status=$?
	printf "return: %d" "$compile_status"
	if (( expect_fail )); then
		if (( compile_status != 0 )); then
			printf " (expected failure)"
		else
			printf " (UNEXPECTED SUCCESS)"
			set -e
			return 1
		fi
		if [[ -e "$out" ]]; then
			echo
			echo "unexpected artifact: $out"
			set -e
			return 1
		fi
	else
		if (( compile_status == 0 )); then
			printf " (expected success)"
		else
			printf " (UNEXPECTED FAILURE)"
			set -e
			return 1
		fi
		if [[ ! -e "$out" ]]; then
			echo
			echo "missing artifact: $out"
			set -e
			return 1
		fi
	fi
	printf "\n\n"
	set -e
}

compare_dags() {
	local left="$1"
	local right="$2"
	COMPARE_COUNT=$((COMPARE_COUNT + 1))
	echo "compare: $left == $right"
	PYTHONPATH=python python - <<'PY' "$left" "$right"
import json, sys
left = json.load(open(sys.argv[1]))
right = json.load(open(sys.argv[2]))
if left != right:
    raise SystemExit(1)
PY
	echo "result: ok"
	echo
}

print_summary() {
	echo "summary:"
	echo "  unit tests: passed"
	echo "  task files checked: $TASK_COUNT"
	echo "  workflow files checked: $WF_COUNT"
	echo "  dag equality checks: $COMPARE_COUNT"
}

if (( $# > 0 )); then
	if [[ $1 =~ \.swl ]]; then
		evaluate_wf $1
	else
		evaluate_task $1
	fi
else
	run_unit_tests
	mkdir -p tests/dag && rm -f tests/dag/*
	for task in tests/*.sh; do
		evaluate_task $task
	done
	for swl in tests/*.swl; do
		evaluate_wf $swl
	done
	compare_dags tests/dag/pipe.json tests/dag/function.json
	compare_dags tests/dag/explicit.json tests/dag/function.json
	mkdir -p tests/cwl
	PYTHONPATH=python python -m swl.transpile.cwl.cli tests/dag/function.json -o tests/cwl/function.pack.cwl
	PYTHONPATH=python python -m swl.transpile.cwl.cli tests/dag/pipe.json -o tests/cwl/pipe.pack.cwl
	PYTHONPATH=python python -m swl.transpile.cwl.cli tests/dag/explicit.json -o tests/cwl/explicit.pack.cwl
	print_summary
fi

