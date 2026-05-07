#!/bin/bash

set -e

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
		swl.ir.test_force_codec
	printf "unit tests passed\n\n"
}

evaluate_task() {
	echo "file: $1"
	echo "syntax:"
	PYTHONPATH=python python -m swl.eval_task $1
	printf "return: $?\n"
	echo "semantic:"
	PYTHONPATH=python python -m swl.eval_task_semantic $1
	printf "return: $?\n\n"
}

evaluate_wf() {
	echo "file: $1"
	echo "syntax:"
	PYTHONPATH=python python -m swl.eval $1
	printf "return: $?\n"
	echo "semantic:"
	PYTHONPATH=python python -m swl.eval_wf_semantic $1
	printf "return: $?\n"
	echo "ir:"
	PYTHONPATH=python python -m swl.eval_ir $1
	printf "return: $?\n"
	echo "force:"
	PYTHONPATH=python python -m swl.eval_force $1
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

