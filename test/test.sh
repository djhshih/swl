#!/bin/bash

evaluate() {
	echo "file: $1"
	PYTHONPATH=../python python -m swl.eval $1
	printf "return: $?\n\n"
}

if (( $# > 0 )); then
	evaluate $1
else
	for swl in *.swl; do
		evaluate $swl
	done
fi

