#!/bin/bash

set -e

run_unit_tests() {
	echo "running unit tests"
	PYTHONPATH=python python -m unittest \
		swl.syntax.lexer \
		swl.syntax.parser \
		swl.syntax.task.test_parser \
		swl.syntax.task.test_interpolation \
		swl.syntax.task.test_bash
	printf "unit tests passed\n\n"
}

evaluate_task() {
	echo "file: $1"
	PYTHONPATH=python python -m swl.eval_task $1
	printf "return: $?\n\n"
}

evaluate_wf() {
	echo "file: $1"
	PYTHONPATH=python python -m swl.eval $1
	printf "return: $?\n\n"
}

if (( $# > 0 )); then
	if [[ $1 =~ \.swl ]]; then
		evaluate_wf $1
	else
		evaluate_task $1
	fi
else
	run_unit_tests
	for task in tests/*.sh; do
		evaluate_task $task
	done
	for swl in tests/*.swl; do
		evaluate_wf $swl
	done
fi

