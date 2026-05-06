# Simple workflow language (SWL)

SWL aims to be easy language for specifying computational workflows by
making use of existing scripts.

We first annotate existing scripts with input and output parameters,
and we can then import them into SWL as callable functions:

```
align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

align | sort | call
```

Here, `|` is the function chain operator that passes records between functions,
and `//` is the record update operator that merges records (right overrides left).
See `doc/README.md` for more details.

This is a work in progress (see `TODO`).

