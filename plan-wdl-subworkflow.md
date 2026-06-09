# Plan: Fix sprocket sub-workflow failures

## Problem

Sprocket 0.26.0 only allows one `workflow` per source file. The WDL emitter currently generates a separate sub-workflow (`main_call_variant`) for mapped/scattered workflow-type steps, plus the main `workflow main`. All three failing workflows (panel, map, map_by) emit two workflows.

## Root Cause

In `emit.py`, `transpile_dag_dict()` (line 76) calls `_subworkflow_to_wdl()` for any DAG step with `type='workflow'`. This generates a standalone `workflow main_call_variant { ... }` definition. Then `_dag_to_wdl` emits the main workflow which calls the sub-workflow inside a scatter block.

## Fix: Inline sub-workflow tasks for mapped steps

When a workflow-type step has a `map` attribute (i.e., it's called inside a scatter), don't generate a separate workflow. Instead, emit the sub-workflow's internal steps (tasks) directly inside the scatter block.

### Example (panel)

**Before:**
```wdl
workflow main_call_variant {
  call align { input: fastq1 = fastq1, ... }
  call sort { input: bam = align.bam, outbase = outbase }
  call call_ { input: bam = sort.bam, ... }
}

workflow main {
  scatter (call_variant_i in range(length(fastq1))) {
    call main_call_variant as call_variant {
      input: fastq1 = fastq1[call_variant_i], ...
    }
  }
  call merge { input: bcf = call_variant.bcf, outbase = "merged" }
  output { File bcf = merge.bcf_out }
}
```

**After:**
```wdl
workflow main {
  scatter (call_variant_i in range(length(fastq1))) {
    call align {
      input: fastq1 = fastq1[call_variant_i], fastq2 = fastq2[call_variant_i], ...
    }
    call sort {
      input: bam = align.bam, outbase = outbase[call_variant_i]
    }
    call call_ {
      input: bam = sort.bam, outbase = outbase[call_variant_i], ref = ref[call_variant_i], ref_fai = ref_fai[call_variant_i]
    }
  }
  call merge { input: bcf = call_.bcf, outbase = "merged" }
  output { File bcf = merge.bcf_out }
}
```

WDL automatically collects variables assigned inside a scatter into arrays, so `call_.bcf` outside the scatter is `Array[File]` — which is what `merge` expects.

## Changes needed

### 1. `emit.py` — `transpile_dag_dict()` (lines 73-80)

Skip sub-workflow generation for workflow-type steps that have a `map` attribute. Instead, register the sub-workflow's internal tasks:

```python
for step in dag.steps:
    if step.id not in tasks:
        if step.type == 'workflow' and getattr(step, 'map', None) is None:
            wdl = _subworkflow_to_wdl(step, workflow_id)
            tasks[step.id] = _wf_name(f'{workflow_id}_{step.id}')
            sub_workflows.append(wdl)
        else:
            tasks[step.id] = _task_name(step.id)
```

For mapped workflow steps, the internal task IDs from `step.definition.dag.steps` need to be registered as tasks. The workflow step itself doesn't get a task entry (it's inlined).

### 2. `emit.py` — `_dag_to_wdl()` (lines 292-314)

When iterating steps, for a mapped workflow-type step, skip the normal call emission (it's handled by `_mapped_step_to_wdl`). Currently the check at line 293-298 already routes mapped steps to `_mapped_step_to_wdl`, so no change needed here.

### 3. `emit.py` — `_mapped_step_to_wdl()` (new logic)

For a workflow-type step, instead of emitting a single call to the sub-workflow, emit calls to each internal step, composing bindings:

- For each internal step, combine its bindings with the sub-workflow's call bindings
- Internal bindings that reference sub-workflow inputs → substitute the scatter expression (e.g., `fastq1[call_variant_i]`)
- Internal bindings that reference other internal steps → keep as-is (same scatter scope)
- Input bindings that are literals → keep as-is

Key insight: the sub-workflow has scalar inputs (`File fastq1`). The outer scatter indexes arrays to produce scalars (`fastq1[call_variant_i]` → `File`). The internal tasks also take scalar inputs. So composing bindings means: wherever an internal task's binding references a sub-workflow input name, emit the scatter-indexed expression instead.

### 4. `emit.py` — Output handling

For mapped workflow-type steps, the sub-workflow's outputs are defined in `step.definition.dag.outputs`. These are referenced by later steps (e.g., `merge` uses `call_variant.bcf`). When inlined, the last internal step's output (e.g., `call_.bcf`) is automatically collected to an array by the scatter. The downstream step binding already references `call_variant.bcf` — but now the step name in the scatter is different. Need to ensure the output name mapping is correct.

Actually, looking at the current panel WDL:
```wdl
call merge { input: bcf = call_variant.bcf, ... }
```

Here `call_variant` is the alias of the sub-workflow call inside the scatter. When inlined, the sub-workflow's internal step `call_` (aliased as `call_` in the WDL) produces `.bcf`. After the scatter, this is referenced as `call_.bcf`.

But the merge step references `call_variant.bcf`. So the output references need to be remapped.

**Approach**: After inlining, track the output mapping: `{workflow_step_id.output_name} → {last_sub_step_id.output_name}`. Then in `_dag_to_wdl`, when emitting the merge step (or any downstream step), substitute `call_variant.bcf` with `call_.bcf`.

This is similar to how `output_renames` works for input/output name collisions, but for the workflow-to-inlined-task mapping.

### 5. Task deduplication

The internal tasks (align, sort, call) may already be defined in the main DAG's task list. The current code at lines 91-94 emits task definitions for all task-type steps. We need to avoid duplicate task definitions.

## Files to modify

- `python/swl/transpile/wdl/emit.py`

## Test updates

- `tests/unit/swl/transpile/wdl/test_emit.py` — update if any assertions on sub-workflow output change
- Integration test: panel, map, map_by should pass with sprocket

## Complexity

Medium-high. The main challenge is binding composition (step 3) and output reference remapping (step 4). The approach is well-defined but requires careful implementation to handle all three workflow types (panel, map, map_by).

## Alternative considered

**Use womtool/cromwell as sprocket fallback**: Rejected since user wants sprocket-only testing. Also womtool 86 doesn't support sub-workflows (WDL 1.1 feature).
