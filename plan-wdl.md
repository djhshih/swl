# Plan: Add WDL Integration Tests

## Goal

Add integration testing for WDL transpilation, analogous to the existing CWL integration tests at `tests/integration/cwl/`. The WDL tests should:

1. Compile SWL workflows to DAG JSON
2. Transpile DAG JSON to WDL 1.1
3. Run each WDL workflow through Cromwell
4. Report PASS/FAIL based on Cromwell exit code

## Background

- WDL transpilation is mature: `python/swl/transpile/wdl/emit.py` (571 lines)
- 15 unit tests in `tests/unit/swl/transpile/wdl/test_emit.py`
- Compile tests generate golden WDL files in `tests/compile/wdl/` (function, pipe, explicit, panel, map, map_by)
- **Cromwell 86** is available at `/home/davids/local/bin/cromwell`
- **womtool 86** is available at `/home/davids/local/bin/womtool`
- No WDL integration tests exist yet
- All input files are stubs (empty/simulated) — same as CWL tests
- Python tasks (`align.sh`, `sort.sh`, `call.sh`, `merge_bcf.sh`) use bash commands (`cp`, `echo`) that work on simulated files

## Implementation Steps

### Step 1: Create directory structure

```
tests/integration/wdl/
├── .gitignore          # ignore *.wdl, outputs/, *.json
├── run.sh              # main test script
├── inputs/             # stub input files (same as CWL: sample1.r1.fq, etc.)
├── jobs/               # JSON input files for Cromwell (one per workflow)
└── outputs/            # Cromwell output dirs (gitignored)
```

### Step 2: Create input stubs

Same input files as CWL integration tests (create in `inputs/`):

```
sample1.r1.fq  sample1.r2.fq  sample2.r1.fq  sample2.r2.fq
ref.fa         ref.fa.amb     ref.fa.ann     ref.fa.bwt
ref.fa.fai     ref.fa.pac     ref.fa.sa      empty.fq
```

Each file contains just the filename as content (same as CWL).

### Step 3: Create WDL job input JSON files

Cromwell accepts JSON input files. For each workflow, create a JSON matching the WDL workflow inputs.

#### `jobs/function.json`
Scalar inputs matching `function.wdl`'s `main` workflow:
```json
{
  "main.fastq1": {"class": "File", "path": "ABS/inputs/sample1.r1.fq"},
  "main.fastq2": {"class": "File", "path": "ABS/inputs/sample1.r2.fq"},
  "main.outbase": "test",
  "main.ref": {"class": "File", "path": "ABS/inputs/ref.fa"},
  "main.ref_amb": {"class": "File", "path": "ABS/inputs/ref.fa.amb"},
  "main.ref_ann": {"class": "File", "path": "ABS/inputs/ref.fa.ann"},
  "main.ref_bwt": {"class": "File", "path": "ABS/inputs/ref.fa.bwt"},
  "main.ref_fai": {"class": "File", "path": "ABS/inputs/ref.fa.fai"},
  "main.ref_pac": {"class": "File", "path": "ABS/inputs/ref.fa.pac"},
  "main.ref_sa": {"class": "File", "path": "ABS/inputs/ref.fa.sa"}
}
```

#### `jobs/panel.json`
Array inputs for scatter + merge. Same structure as function inputs but with arrays for scattered columns.
```json
{
  "main.fastq1": [{"class": "File", "path": "..."}, {"class": "File", "path": "..."}],
  "main.fastq2": [...],
  ...
}
```

#### `jobs/map.json`
WDL `map.wdl` has workflow `main` with input `Array[Call_variantInput] xs` (a struct array). Need to create an array of structs.
```json
{
  "main.xs": [
    {
      "fastq1": {"class": "File", "path": "..."},
      "fastq2": {"class": "File", "path": "..."},
      ...
    },
    {
      "fastq1": {"class": "File", "path": "..."},
      ...
    }
  ]
}
```

> **Question**: Cromwell uses its own File class `{"class": "File", "path": "..."}`. The WDL compile test generates `function.wdl` which uses `String` and `File` types directly. The JSON must use Cromwell's schema, NOT CWL's. Need to verify what Cromwell accepts.

#### `jobs/map_by.json`
WDL `map_by.wdl` has workflow `main` with `Array[String] outbase, Array[File] fastq1, ...` input. Provide arrays of values:
```json
{
  "main.outbase": ["sample1", "sample2"],
  "main.fastq1": [{"class": "File", "path": "..."}, {"class": "File", "path": "..."}],
  ...
}
```

### Step 4: Write `run.sh`

The test script should:

1. Create stub input files
2. Compile each SWL to DAG JSON (`swl.compile`)
3. Transpile each DAG to WDL (`swl.transpile.wdl`)
4. For each workflow, run with Cromwell:
   ```bash
   cromwell run wdl_file -i inputs.json
   ```
5. Extract exit code and report PASS/FAIL

**Important considerations:**

- Cromwell is a Java jar. Running `cromwell run` may be slow (JVM startup). Use `--type WDL --type-version 1.1` if needed.
- Cromwell writes outputs to a timestamped directory. The test script should capture the output JSON and check for success.
- Cromwell returns exit code 0 on success, non-zero on failure.
- The JS `inputs_file` paths need to be absolute. The script should resolve them.
- Test with `function` first (simplest, scalar inputs), then `panel`, `map`, `map_by`.

**Test selection:**
Start with the same workflows as CWL:
- `function` (scalar, simple pipe)
- `panel` (map + merge)
- `map` (scatter over struct array)
- `map_by` (grouped scatter)

**Note on Cromwell exit code:**
Cromwell prints "WORKFLOW FINISHED WITH STATUS: Succeeded" on success. The exit code should reflect success. If exit code alone is insufficient, parse stdout for the status string.

### Step 5: `womtool validate` as a pre-check

Before running Cromwell, add a `womtool validate` step:

```bash
womtool validate wdl_file
```

This gives a fast syntax/validation check before the slow Cromwell execution. Include this in the test flow but only count Cromwell execution for PASS/FAIL.

### Step 6: Integrate with master `test.sh`

Update `/home/davids/projects/swl/test.sh` to add:

```bash
run_suite "integration-wdl" "tests/integration/wdl/run.sh"
```

### Step 7: Create `.gitignore`

```
outputs/
*.wdl
*.json
```

### Step 8: Add SWL symlinks

If needed, symlink SWL files from `tests/integration/swl/` (which is shared with CWL). The CWL integration tests already have SWL files at `tests/integration/swl/`. The WDL test can share these.

### Step 9: Update `tests/compile/test.sh`

The compile tests currently transpile WDL golden files but don't test `map_by` for WDL compilation (since map_by WDL support was added later). Verify the compile test already includes `map_by.wdl`:

Check the compile test at `tests/compile/test.sh` lines 143-153. If `map_by` is not in the WDL transpile step, add it. Looking at the compile test — it transpiles `map_by.json` to WDL on line 149 (`../wdl/map_by.wdl`). The file comparison is `pipe.wdl == function.wdl` and `explicit.wdl == function.wdl`. This doesn't test `map_by.wdl` content, just that it transpiles without error. **No change needed.**

## Risks and Open Questions

1. **Cromwell startup time**: JVM startup is slow, adding ~5-10s per test. With 4 workflows, that's 20-40s. Acceptable for integration tests.

2. **Test input files are empty/stubs**: The actual bash commands in the tasks (`bwa mem`, `samtools`, etc.) will FAIL on empty inputs. The CWL tests use `cp` commands via shell interpolation that work on empty files. Need to verify the WDL commands (which use the same task scripts) will also work.

   **Check**: What do the task scripts actually contain? The SWL compiles to DAG with `script` or `body` fields containing bash commands like `bwa mem -t ${cpu} ...`. These commands reference the inputs. If inputs exist as empty files, the commands may fail (e.g., `bwa mem` on empty FASTQ).

   **Resolution**: The CWL tests work because the test scripts use `cp` instead of real tools (the SWL shell tasks were modified for testing). Verify the WDL integration tests use the same test-specific task scripts (`tests/compile/swl/align.sh`, etc.) which use `cp` for test purposes.

3. **Cromwell File path format**: Cromwell uses `{"class": "File", "path": "..."}` in JSON inputs. Need absolute paths or relative-to-workflow-root paths.

4. **WDL version**: The transpiler emits `version 1.1`. Cromwell 86 supports WDL 1.1.

5. **Docker containers**: The WDL tasks declare `container` requirements. Cromwell will try to pull these containers. If containers aren't available, tasks will fail. The test scripts should use tasks without container requirements, or ensure containers are pre-pulled.

   **Current tasks check**: Looking at the task scripts in `tests/compile/swl/`:
   - `align.sh`: uses `bwa` and `samtools` (needs containers)
   - `sort.sh`: uses `samtools` (needs containers)
   - `call.sh`: uses `bcftools` (needs containers)
   
   The compiled WDL includes `container` requirements. Cromwell will try to pull these. If containers aren't available, the integration test will fail.

   **Options**:
   a. Use `--no-container` or disable containers in Cromwell options
   b. Pre-pull containers before running tests
   c. Use different task scripts without container requirements

   For initial implementation, use Cromwell's `--options` file to disable container execution.

## Implementation Order

1. Create `tests/integration/wdl/` directory structure and `.gitignore`
2. Create stub input files in `tests/integration/wdl/inputs/`
3. Write `run.sh` with `womtool validate` + `cromwell run` for each workflow
4. Create JSON input files in `tests/integration/wdl/jobs/`
5. Test manually with a single workflow (function)
6. Add remaining workflows (panel, map, map_by)
7. Update master `test.sh`
8. Verify full test suite passes

## Expected Output

After implementation:
```
===== unit =====
238 tests pass
===== compile =====
10 files compile
===== integration =====
CWL: 4/4 pass
WDL: 4/4 pass
===== RESULTS =====
Passed: 4 / 4
```
