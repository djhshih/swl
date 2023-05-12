#!/bin/bash

evaluate() {
	PYTHONPATH=../python python -m swl.eval $1
}

for swl in *.swl; do
	echo "file: $swl"
	evaluate $swl
	printf "return: $?\n\n"
done

