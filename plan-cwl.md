# Plan: Restructure tests/ and add cwltool integration testing

## Goal

1. Restructure `tests/` directory: move source files into `tests/compile/`, organize by test category
2. Split `test.sh` into a master script that delegates to sub-runners
3. Add integration tests that run `cwltool` on transpiled CWL using stub task scripts

---

## Background: integration test approach

The integration test `swl/` directory contains symlinks to the real SWL files plus committed stub `.sh` files. Since `import "align.sh"` resolves relative to the `.swl` file's directory, and `os.path.abspath` preserves symlinks, the stubs shadow the originals at compile time. The runner compiles directly from `swl/` — no file copying needed.

---

## Current state

```
swl/
├── test.sh                           # monolithic test runner
└── tests/
    ├── __init__.py
    ├── .gitignore                    # ignores dag/
    ├── align.sh, sort.sh, call.sh, ...   # 7 task scripts
    ├── function.swl, pipe.swl, ...       # 9 SWL files
    ├── dag/                          # compiled DAG JSONs (generated)
    ├── cwl/  (.gitignore: *.cwl)     # transpiled CWL (generated)
    ├── wdl/  (.gitignore: *.wdl)     # transpiled WDL (generated)
    ├── nf/   (.gitignore: *.nf)      # transpiled Nextflow (generated)
    ├── smk/  (.gitignore: *.smk)     # transpiled Snakemake (generated)
    └── unit/                         # Python unit tests (unchanged)
```

---

## Target directory structure

```
swl/
├── test.sh                           # MASTER — calls sub-test.sh scripts
└── tests/
    ├── __init__.py
    ├── compile/                      # source originals + golden tests
    │   ├── .gitignore                # ignores dag/ cwl/ wdl/ nf/ smk/
    │   ├── swl/                      # source files
    │   │   ├── align.sh
    │   │   ├── sort.sh
    │   │   ├── call.sh
    │   │   ├── cat.sh
    │   │   ├── merge_bcf.sh
    │   │   ├── realign.sh
    │   │   ├── bad_outbase.sh
    │   │   ├── function.swl
    │   │   ├── pipe.swl
    │   │   ├── explicit.swl
    │   │   ├── panel.swl
    │   │   ├── map.swl
    │   │   ├── map_by.swl
    │   │   ├── bad_explicit.swl
    │   │   ├── bad_pipe.swl
    │   │   ├── partial.swl
    │   │   └── import_partial.swl
    │   ├── dag/                      # compiled DAG JSONs (generated)
    │   ├── cwl/                      # transpiled CWL (generated)
    │   ├── wdl/                      # transpiled WDL (generated)
    │   ├── nf/                       # transpiled Nextflow (generated)
    │   ├── smk/                      # transpiled Snakemake (generated)
    │   └── test.sh                   # compile + golden comparison tests
    │
    ├── integration/
    │   ├── swl/                      # SWL symlinks + committed stub scripts
    │   │   ├── align.sh              # stub (shared by all language backends)
    │   │   ├── sort.sh               # stub
    │   │   ├── call.sh               # stub
    │   │   ├── function.swl ─────────┐
    │   │   ├── pipe.swl             ├──> ../../compile/swl/<same>
    │   │   ├── explicit.swl         │
    │   │   ├── panel.swl            │
    │   │   └── map.swl ─────────────┘
    │   └── cwl/
    │       ├── .gitignore            # ignores outputs/
    │       ├── function.cwl          # transpiled with stubs (generated)
    │       ├── pipe.cwl              # (generated)
    │       ├── explicit.cwl          # (generated)
    │       ├── panel.cwl             # (generated)
    │       ├── map.cwl               # (generated)
    │       ├── jobs/                 # cwltool job YAML files
    │       │   ├── function.yml
    │       │   ├── panel.yml
    │       │   └── map.yml
    │       ├── inputs/               # mock input files (created by run.sh)
    │       └── outputs/              # cwltool output dir (generated)
    │       └── run.sh                # integration test runner
    │
    └── unit/                         # Python unit tests (unchanged)
        ├── __init__.py
        ├── test.sh                   # unit test runner
        └── swl/                      # (unchanged)

```

---

## Master `test.sh` (at project root)

```bash
#!/bin/bash
set -e

PASS=0
FAIL=0

run_suite() {
    local name="$1" path="$2"
    echo "===== $name ====="
    if bash "$path"; then
        echo "PASS: $name"
        PASS=$((PASS + 1))
    else
        echo "FAIL: $name"
        FAIL=$((FAIL + 1))
    fi
    echo
}

run_suite "unit"       "tests/unit/test.sh"
run_suite "compile"    "tests/compile/test.sh"
run_suite "integration" "tests/integration/cwl/run.sh"

echo "===== RESULTS ====="
echo "Passed: $PASS / $((PASS + FAIL))"
```

---

## `tests/compile/test.sh` (extracted from current `test.sh`)

Ported from the current monolithic `test.sh`. Runs from within `tests/compile/swl/`, outputs to sibling dirs.

```bash
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
```

---

## `tests/unit/test.sh`

Simple wrapper for existing unit test discovery:

```bash
#!/bin/bash
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
echo "running unit tests"
PYTHONPATH="$ROOT/python:$ROOT/tests/unit" python -m unittest \
    tests.unit.swl.syntax.wf.test_lexer \
    tests.unit.swl.syntax.wf.test_parser \
    tests.unit.swl.syntax.task.test_parser \
    tests.unit.swl.syntax.task.test_interpolation \
    tests.unit.swl.syntax.task.test_bash \
    tests.unit.swl.semantic.task.test_type \
    tests.unit.swl.semantic.wf.test_check \
    tests.unit.swl.ir.test_lower \
    tests.unit.swl.ir.test_force \
    tests.unit.swl.ir.test_force_codec \
    tests.unit.swl.transpile.cwl.test_emit \
    tests.unit.swl.transpile.wdl.test_emit \
    tests.unit.swl.transpile.nf.test_emit
printf "unit tests passed\n\n"
```

---

## Integration test runner: `tests/integration/cwl/run.sh`

The runner compiles from `tests/integration/swl/` (SWL symlinks + stub symlinks), then transpiles and runs cwltool.

```bash
#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
INT_DIR="$ROOT/tests/integration/cwl"
SWL_DIR="$ROOT/tests/integration/swl"
JOBS_DIR="$INT_DIR/jobs"
INPUTS_DIR="$INT_DIR/inputs"
OUTPUTS_DIR="$INT_DIR/outputs"
PASS=0
FAIL=0

cwltool --version

# 1. Create mock input files
mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR"
for f in sample1.r1.fq sample1.r2.fq sample2.r1.fq sample2.r2.fq \
         ref.fa ref.fa.amb ref.fa.ann ref.fa.bwt \
         ref.fa.fai ref.fa.pac ref.fa.sa empty.fq; do
  echo "$f" > "$INPUTS_DIR/$f"
done

# 2. Compile SWL → DAG (stub symlinks in swl/ shadow originals at compile time)
PYTHONPATH="$ROOT/python"
for swl in function pipe explicit panel map; do
  python -m swl.compile "$SWL_DIR/${swl}.swl" -o "$INT_DIR/${swl}.json"
done

# 3. Transpile DAG → CWL (stub bodies embedded in the DAG JSONs)
for swl in function pipe explicit panel map; do
  python -m swl.transpile.cwl "$INT_DIR/${swl}.json" -o "$INT_DIR/${swl}.cwl"
done

# 4. Run each workflow through cwltool
run_test() {
  local name="$1"
  local cwl="$INT_DIR/${name}.cwl"
  local job="$JOBS_DIR/${name}.yml"
  local outdir="$OUTPUTS_DIR/$name"
  mkdir -p "$outdir"
  echo "=== $name ==="
  if cwltool --outdir "$outdir" "$cwl" "$job" 2>"$outdir/stderr.log"; then
    echo "PASS: $name"
    PASS=$((PASS + 1))
  else
    echo "FAIL: $name (see $outdir/stderr.log)"
    FAIL=$((FAIL + 1))
  fi
}

run_test function
run_test panel
run_test map

echo "---"
echo "Passed: $PASS / $((PASS + FAIL))"
```

---

## Stub scripts

Committed directly in `tests/integration/swl/` alongside the SWL symlinks.

### `swl/align.sh`
```bash
cp "${fastq1}" "${outbase}.bam"
```

### `swl/sort.sh`
```bash
cp "${bam}" "${outbase}.bam"
cp "${bam}" "${outbase}.bai"
```

### `swl/call.sh`
```bash
cp "${bam}" "${outbase}.bcf"
```

Stubs use `cp` (not `touch`) so cwltool's dependency wiring is exercised — missing inputs cause non-zero exit.

---

## Job YAML files

cwltool resolves paths relative to the job file's directory. Jobs live in `tests/integration/cwl/jobs/`, inputs in `tests/integration/cwl/inputs/`, so paths use `../inputs/<name>`.

### `tests/integration/cwl/jobs/function.yml`
```yaml
fastq1:
  class: File
  path: ../inputs/sample1.r1.fq
fastq2:
  class: File
  path: ../inputs/sample1.r2.fq
outbase: sample1
ref:
  class: File
  path: ../inputs/ref.fa
ref_amb:    { class: File, path: ../inputs/ref.fa.amb }
ref_ann:    { class: File, path: ../inputs/ref.fa.ann }
ref_bwt:    { class: File, path: ../inputs/ref.fa.bwt }
ref_fai:    { class: File, path: ../inputs/ref.fa.fai }
ref_pac:    { class: File, path: ../inputs/ref.fa.pac }
ref_sa:     { class: File, path: ../inputs/ref.fa.sa }
```

### `tests/integration/cwl/jobs/panel.yml`
```yaml
fastq1:
  class: File
  path: ../inputs/sample1.r1.fq
fastq2:
  class: File
  path: ../inputs/sample1.r2.fq
outbase: panel-test
```

### `tests/integration/cwl/jobs/map.yml`
```yaml
xs:
  - fastq1:  { class: File, path: ../inputs/sample1.r1.fq }
    fastq2:  { class: File, path: ../inputs/sample1.r2.fq }
    outbase: sample1
    ref:     { class: File, path: ../inputs/ref.fa }
    ref_amb: { class: File, path: ../inputs/ref.fa.amb }
    ref_ann: { class: File, path: ../inputs/ref.fa.ann }
    ref_bwt: { class: File, path: ../inputs/ref.fa.bwt }
    ref_fai: { class: File, path: ../inputs/ref.fa.fai }
    ref_pac: { class: File, path: ../inputs/ref.fa.pac }
    ref_sa:  { class: File, path: ../inputs/ref.fa.sa }
  - fastq1:  { class: File, path: ../inputs/sample2.r1.fq }
    fastq2:  { class: File, path: ../inputs/sample2.r2.fq }
    outbase: sample2
    ref:     { class: File, path: ../inputs/ref.fa }
    ref_amb: { class: File, path: ../inputs/ref.fa.amb }
    ref_ann: { class: File, path: ../inputs/ref.fa.ann }
    ref_bwt: { class: File, path: ../inputs/ref.fa.bwt }
    ref_fai: { class: File, path: ../inputs/ref.fa.fai }
    ref_pac: { class: File, path: ../inputs/ref.fa.pac }
    ref_sa:  { class: File, path: ../inputs/ref.fa.sa }
```

---

## Validation matrix

| Workflow | Expected outputs | Verification |
|----------|-----------------|--------------|
| `function` | `sample1.bam`, `sample1.bai`, `sample1.bcf` | Files exist; `sample1.bam` content == `sample1.r1.fq` |
| `panel` | `panel-test.bam` | File exists; content == `sample1.r1.fq` |
| `map` | `sample1.{bam,bai,bcf}`, `sample2.{bam,bai,bcf}` | All 6 exist; content matches respective fastq1 |

---

## Implementation steps

1. Create `tests/compile/swl/` and `tests/compile/{dag,cwl,wdl,nf,smk}/`
2. Move all `tests/*.swl` → `tests/compile/swl/` (9 files)
3. Move all `tests/*.sh` → `tests/compile/swl/` (7 files)
4. Move `tests/dag/` → `tests/compile/dag/`
5. Move `tests/cwl/` → `tests/compile/cwl/`
6. Move `tests/wdl/` → `tests/compile/wdl/`
7. Move `tests/nf/` → `tests/compile/nf/`
8. Move `tests/smk/` → `tests/compile/smk/`
9. Write `tests/compile/.gitignore` (ignore `dag/ cwl/ wdl/ nf/ smk/`)
10. Write `tests/compile/test.sh`
11. Write `tests/unit/test.sh`
12. Create `tests/integration/swl/` with SWL symlinks → `compile/swl/` + committed stub scripts
13. Create `tests/integration/cwl/` directory structure
14. Write job YAML files in `tests/integration/cwl/jobs/`
15. Write `tests/integration/cwl/.gitignore` (ignore `outputs/`)
16. Write `tests/integration/cwl/run.sh`
17. Write master `test.sh` at project root
18. Run master `test.sh` to verify all suites pass
