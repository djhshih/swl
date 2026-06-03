import os
import tempfile
import unittest as ut

from swl.semantic.wf import type as wf_type
from swl.semantic.wf.check import Checker
from swl.semantic.wf.infer import ClosedRecord, ClosureValue, ComputationValue, FunctionValue, UnknownValue, application_result, apply_function


_ALIGN = '''# @ Align
# in
#   fastq1, fastq2 file
#   ref file
#   ref_fai file
#   outbase str
# out
#   bam file = ${outbase}.bam
# run
#   cpu = 2
echo align
'''

_SORT = '''# @ Sort
# in
#   bam file
#   outbase str
# out
#   bam file = ${outbase}.bam
#   bai file = ${outbase}.bai
echo sort
'''

_CALL = '''# @ Call
# in
#   bam file
#   ref file
#   ref_fai file
#   outbase str
# out
#   bcf file = ${outbase}.bcf
echo call
'''

_BAD_OUTBASE = '''# @ Bad
# out
#   outbase file = result.txt
echo bad
'''

_PIPE = '''align = import "align.sh"
sort = import "sort.sh"
call = import "call.sh"
align | sort | call
'''

_FUNCTION = '''align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

\\x ->
    a = align x
    s = sort ( x // a )
    c = call ( x // a // s )
    a // s // c
'''

_EXPLICIT = '''align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

\\x ->
    a = align x
    s = sort { bam: a.bam, outbase: x.outbase }
    c = call { bam: s.bam, ref: x.ref, ref_fai: x.ref_fai, outbase: x.outbase }
    { bam: s.bam, bai: s.bai, bcf: c.bcf }
'''

_PARTIAL = '''align = import "align.sh"
align_hg38 = align {
  ref: "hg38.fa",
  ref_fai: "hg38.fa.fai",
  cpu: 2
}
align_hg38
'''

_IMPORT_PARTIAL = '''partial = import "partial.swl"
partial
'''

_LAMBDA_APPLY = '''mk = \\x -> { y: x.a }
\\z ->
    r = mk z
    { y: r.y }
'''

_RECORD_WORKFLOW = '''\\x ->
    { foo: x.foo, bar: x.bar }
'''

_IMPORT_RECORD = '''recorder = import "record.swl"
recorder
'''

_TASK_RESULT_WORKFLOW = '''align = import "align.sh"
\\x ->
    align x
'''

_IMPORT_TASK_RESULT = '''w = import "task_result.swl"
w
'''

_INPUT_PROP_WORKFLOW = '''align = import "align.sh"
\\x ->
    a = align x
    { bam: a.bam }
'''

_IMPORT_INPUT_PROP = '''w = import "input_prop.swl"
w
'''

_PARTIAL_INNER_WORKFLOW = '''align = import "align.sh"
\\x ->
    f = align { ref: "hg38.fa", ref_fai: "hg38.fa.fai" }
    f x
'''

_IMPORT_PARTIAL_INNER = '''w = import "partial_inner.swl"
w
'''

_BATCH_WORKFLOW_PARTIAL = '''mkp = import "mk_align_partial.swl"
merge = import "merge.sh"
\\xs ->
    ys = map mkp xs
    merge { bam: ys.bam, outbase: "merged" }
'''

_CYCLE_A = '''b = import "b.swl"
b
'''

_CYCLE_B = '''a = import "a.swl"
a
'''

_CHAINABLE = '''# @ Chainable
# in
#   bam file
# out
#   bcf file = out.bcf
echo chainable
'''

_WORKFLOW_CHAIN = '''w = import "task_result.swl"
call = import "chainable.sh"
w | call
'''

_WORKFLOW_APPLY = '''w = import "task_result.swl"
\\x ->
    y = w x
    { bam: y.bam }
'''

_SCALAR_APPLY = '''align = import "align.sh"
\\x ->
    result = align "reads_1.fq"
    { bam: result.bam }
'''

_BAD_PIPE = '''align = import "bad_outbase.sh"
sort = import "sort.sh"
align | sort
'''

_BAD_EXPLICIT_PARTIAL_FIELD = '''align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

\\x ->
    a = align x
    s = sort { bam: a.bam, outbase: x.outbase }
    c = call { bam: s.bam, ref: x.ref, outbase: x.outbase }
    { bcf: c.bcf }
'''

_DUP_BIND = '''x = 1
x = 2
\\y -> y
'''

_DUP_LAMBDA_PARAM = '''x = 1
\\x -> x
'''

_DUP_IMPORT = '''x = import "align.sh"
x = import "sort.sh"
\\y -> y
'''

_LAMBDA_BODY_SHADOW = '''x = 1
\\y ->
    x = 2
    x
'''

_NON_FUNCTION_RECORD = '''{ foo: 1 }
'''

_FORWARD_REF = '''inner = import "align.sh"
x = inner y
y = { fastq1: "a", fastq2: "b", ref: "r", ref_fai: "r.fai", outbase: "o" }
\\z -> z
'''

_SELF_REF = '''inner = import "align.sh"
x = inner x
\\z -> z
'''

_VALID_CROSS_REF = '''inner = import "align.sh"
x = { fastq1: "a", fastq2: "b", ref: "r", ref_fai: "r.fai", outbase: "o" }
y = inner x
\\z -> z
'''

_NON_FUNCTION_SCALAR = '''1
'''

_NON_FUNCTION_APPLY = '''align = import "align.sh"
align {
  fastq1: "r1.fq",
  fastq2: "r2.fq",
  ref: "hg38.fa",
  ref_fai: "hg38.fa.fai",
  outbase: "sample"
}
'''

_MERGE = '''# @ Merge
# in
#   bam [file]
#   outbase str
# out
#   bam file = ${outbase}.bam
echo merge
'''

_BATCH_OK = '''align = import "align.sh"
merge = import "merge.sh"
\\xs ->
    calls = map align xs
    merge { bam: calls.bam, outbase: "merged" }
'''

_BATCH_BAD_FIELD = '''align = import "align.sh"
\\xs ->
    calls = map align xs
    { x: calls.nope }
'''

_BATCH_NON_ARRAY = '''align = import "align.sh"
\\x ->
    calls = map align x
    { bam: calls.bam }
'''

_BATCH_INNER = '''align = import "align.sh"
\\xs ->
    calls = map align xs
    { bam: calls.bam }
'''

_BATCH_ON_BATCH = '''inner = import "batch_inner.swl"
\\xs ->
    ys = map inner xs
    { bam: ys.bam }
'''

_BATCH_NON_FUNCTION = '''\\xs ->
    ys = map { foo: xs } xs
    { y: ys.foo }
'''


class TestWorkflowCheck(ut.TestCase):
    def _write(self, root, name, content):
        path = os.path.join(root, name)
        with open(path, 'w') as f:
            f.write(content)
        return path

    def _make_fixture_dir(self):
        td = tempfile.TemporaryDirectory()
        self.addCleanup(td.cleanup)
        root = td.name
        self._write(root, 'align.sh', _ALIGN)
        self._write(root, 'sort.sh', _SORT)
        self._write(root, 'call.sh', _CALL)
        self._write(root, 'bad_outbase.sh', _BAD_OUTBASE)
        self._write(root, 'merge.sh', _MERGE)
        return root

    def test_load_pipe_workflow(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'pipe.swl', _PIPE)
        result = Checker().load(path)
        self.assertEqual(sorted(result.imports.keys()), ['align', 'call', 'sort'])
        self.assertEqual(result.errors, [])
        self.assertEqual(result.inferred_inputs, set())
        self.assertIn('bcf', result.signature.outputs)

    def test_load_function_workflow_infers_inputs(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'function.swl', _FUNCTION)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertIn('fastq1', result.inferred_inputs)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)
        self.assertIn('bcf', result.signature.outputs)

    def test_load_explicit_workflow_infers_inputs(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'explicit.swl', _EXPLICIT)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertIn('fastq1', result.inferred_inputs)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)
        self.assertIn('bcf', result.signature.outputs)

    def test_import_workflow_signature(self):
        root = self._make_fixture_dir()
        self._write(root, 'partial.swl', _PARTIAL)
        path = self._write(root, 'import_partial.swl', _IMPORT_PARTIAL)
        result = Checker().load(path)
        self.assertIn('partial', result.imports)
        self.assertEqual(result.imports['partial'].kind, 'workflow')
        self.assertIn('bam', result.imports['partial'].signature.outputs)

    def test_partial_application_returns_a_function_signature(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'partial.swl', _PARTIAL)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertIn('fastq1', result.signature.inputs)
        self.assertIn('fastq2', result.signature.inputs)
        self.assertIn('outbase', result.signature.inputs)
        self.assertNotIn('ref', result.signature.inputs)
        self.assertNotIn('ref_fai', result.signature.inputs)
        self.assertIn('bam', result.signature.outputs)

    def test_closure_preserves_bound_value_fields(self):
        root = self._make_fixture_dir()
        checker = Checker()
        imported = checker._load_import('align', os.path.join(root, 'align.sh'))
        fn = FunctionValue('align', imported.signature, 'task')
        closure = apply_function(checker, fn, ClosedRecord({'ref': object(), 'ref_fai': object()}), set())
        self.assertIsInstance(closure, ClosureValue)
        self.assertIn('ref', closure.bound_value.fields)
        self.assertIn('ref_fai', closure.bound_value.fields)
        self.assertNotIn('ref', closure.signature.inputs)
        self.assertNotIn('ref_fai', closure.signature.inputs)

    def test_lambda_application_evaluates_body_for_signature(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'lambda_apply.swl', _LAMBDA_APPLY)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertEqual(sorted(result.signature.inputs.keys()), ['a'])
        self.assertEqual(sorted(result.signature.outputs.keys()), ['y'])

    def test_computation_preserves_output_record_shape(self):
        root = self._make_fixture_dir()
        checker = Checker()
        imported = checker._load_import('align', os.path.join(root, 'align.sh'))
        fn = FunctionValue('align', imported.signature, 'task')
        value = application_result(
            checker, fn,
            ClosedRecord({
                'fastq1': UnknownValue(),
                'fastq2': UnknownValue(),
                'ref': UnknownValue(),
                'ref_fai': UnknownValue(),
                'outbase': UnknownValue(),
            }),
            set(),
            [],
        )
        self.assertIsInstance(value, ComputationValue)
        self.assertIn('bam', value.output_value.fields)

    def test_import_workflow_record_signature(self):
        root = self._make_fixture_dir()
        self._write(root, 'record.swl', _RECORD_WORKFLOW)
        path = self._write(root, 'import_record.swl', _IMPORT_RECORD)
        result = Checker().load(path)
        self.assertIn('recorder', result.imports)
        self.assertIn('foo', result.imports['recorder'].signature.inputs)
        self.assertIn('bar', result.imports['recorder'].signature.inputs)
        self.assertIn('foo', result.imports['recorder'].signature.outputs)
        self.assertIn('bar', result.imports['recorder'].signature.outputs)

    def test_import_workflow_task_result_signature(self):
        root = self._make_fixture_dir()
        self._write(root, 'task_result.swl', _TASK_RESULT_WORKFLOW)
        path = self._write(root, 'import_task_result.swl', _IMPORT_TASK_RESULT)
        result = Checker().load(path)
        self.assertIn('w', result.imports)
        self.assertIn('bam', result.imports['w'].signature.outputs)
        self.assertIn('fastq1', result.imports['w'].signature.inputs)
        self.assertIn('outbase', result.imports['w'].signature.inputs)

    def test_import_workflow_input_propagation(self):
        root = self._make_fixture_dir()
        self._write(root, 'input_prop.swl', _INPUT_PROP_WORKFLOW)
        path = self._write(root, 'import_input_prop.swl', _IMPORT_INPUT_PROP)
        result = Checker().load(path)
        self.assertIn('w', result.imports)
        self.assertIn('fastq1', result.imports['w'].signature.inputs)
        self.assertIn('fastq2', result.imports['w'].signature.inputs)
        self.assertIn('ref', result.imports['w'].signature.inputs)
        self.assertIn('ref_fai', result.imports['w'].signature.inputs)
        self.assertIn('outbase', result.imports['w'].signature.inputs)

    def test_imported_workflow_with_inner_partial_preserves_remaining_inputs_only(self):
        root = self._make_fixture_dir()
        self._write(root, 'partial_inner.swl', _PARTIAL_INNER_WORKFLOW)
        path = self._write(root, 'import_partial_inner.swl', _IMPORT_PARTIAL_INNER)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertIn('w', result.imports)
        self.assertEqual(sorted(result.imports['w'].signature.inputs.keys()), ['fastq1', 'fastq2', 'outbase'])
        self.assertEqual(
            {k: getattr(getattr(v, 'type', None), 'value', None) for k, v in result.imports['w'].signature.inputs.items()},
            {'fastq1': 'file', 'fastq2': 'file', 'outbase': 'str'},
        )

    def test_batch_imported_workflow_with_inner_partial_propagates_concrete_table_columns(self):
        root = self._make_fixture_dir()
        self._write(root, 'mk_align_partial.swl', _PARTIAL_INNER_WORKFLOW)
        path = self._write(root, 'batch_workflow_partial.swl', _BATCH_WORKFLOW_PARTIAL)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertTrue(isinstance(result.root_input_type, wf_type.TableType))
        self.assertEqual(
            result.root_input_type.columns,
            {
                'fastq1': wf_type.FILE,
                'fastq2': wf_type.FILE,
                'outbase': wf_type.STR,
            },
        )
        self.assertEqual(sorted(result.signature.inputs.keys()), ['fastq1', 'fastq2', 'outbase'])

    def test_circular_workflow_import_fails(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'a.swl', _CYCLE_A)
        self._write(root, 'b.swl', _CYCLE_B)
        with self.assertRaises(ValueError):
            Checker().load(path)

    def test_direct_workflow_import_in_chain(self):
        root = self._make_fixture_dir()
        self._write(root, 'chainable.sh', _CHAINABLE)
        self._write(root, 'task_result.swl', _TASK_RESULT_WORKFLOW)
        path = self._write(root, 'workflow_chain.swl', _WORKFLOW_CHAIN)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertIn('call', result.imports)
        self.assertEqual(result.imports['w'].kind, 'workflow')

    def test_direct_workflow_import_in_application(self):
        root = self._make_fixture_dir()
        self._write(root, 'task_result.swl', _TASK_RESULT_WORKFLOW)
        path = self._write(root, 'workflow_apply.swl', _WORKFLOW_APPLY)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertIn('fastq1', result.inferred_inputs)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)
        self.assertIn('bam', result.signature.outputs)

    def test_scalar_task_application_lifts_to_first_input_but_field_access_on_partial_is_an_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'scalar_apply.swl', _SCALAR_APPLY)
        result = Checker().load(path)
        self.assertIsNotNone(result.signature)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)
        self.assertNotIn('fastq1', result.inferred_inputs)
        self.assertTrue(any('Cannot access field on function value' in err for err in result.errors))

    def test_bad_pipe_reports_chain_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'bad_pipe.swl', _BAD_PIPE)
        result = Checker().load(path)
        self.assertEqual(len(result.errors), 1)
        self.assertIn('outbase', result.errors[0])

    def test_field_access_on_partial_function_reports_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'bad_explicit_partial_field.swl', _BAD_EXPLICIT_PARTIAL_FIELD)
        result = Checker().load(path)
        self.assertTrue(any('Cannot access field on function value' in err for err in result.errors))

    def test_workflow_record_final_value_reports_not_a_function(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'record_final.swl', _NON_FUNCTION_RECORD)
        result = Checker().load(path)
        self.assertIsNone(result.signature)
        self.assertIn('Workflow must evaluate to a function', result.errors)

    def test_workflow_scalar_final_value_reports_not_a_function(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'scalar_final.swl', _NON_FUNCTION_SCALAR)
        result = Checker().load(path)
        self.assertIsNone(result.signature)
        self.assertIn('Workflow must evaluate to a function', result.errors)

    def test_duplicate_binding_in_scope_reports_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'dup_bind.swl', _DUP_BIND)
        result = Checker().load(path)
        self.assertIn('Duplicate binding in scope: x', result.errors)

    def test_lambda_param_may_shadow_outer_scope(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'dup_lambda.swl', _DUP_LAMBDA_PARAM)
        result = Checker().load(path)
        self.assertNotIn('Duplicate binding in scope: x', result.errors)

    def test_duplicate_import_binding_reports_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'dup_import.swl', _DUP_IMPORT)
        result = Checker().load(path)
        self.assertIn('Duplicate binding in scope: x', result.errors)

    def test_lambda_body_may_shadow_outer_scope(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'lambda_body_shadow.swl', _LAMBDA_BODY_SHADOW)
        result = Checker().load(path)
        self.assertNotIn('Duplicate binding in scope: x', result.errors)

    # P3: forward-reference detection ----------------------------------------

    def test_forward_reference_reports_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'fwd_ref.swl', _FORWARD_REF)
        result = Checker().load(path)
        self.assertTrue(
            any('Forward reference' in err for err in result.errors),
            f'Expected forward reference error, got: {result.errors}',
        )

    def test_self_reference_reports_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'self_ref.swl', _SELF_REF)
        result = Checker().load(path)
        self.assertTrue(
            any('references itself' in err for err in result.errors),
            f'Expected self-reference error, got: {result.errors}',
        )

    def test_prior_binding_reference_is_valid(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'valid_cross.swl', _VALID_CROSS_REF)
        result = Checker().load(path)
        self.assertEqual([], result.errors)

    def test_batch_map_workflow_infers_concrete_table_inputs(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'batch_ok.swl', _BATCH_OK)
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertEqual(sorted(result.signature.inputs.keys()), ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'])
        self.assertIn('bam', result.signature.outputs)
        self.assertTrue(isinstance(result.workflow_type, wf_type.FunctionType))
        self.assertTrue(isinstance(result.root_input_type, wf_type.TableType))
        self.assertEqual(
            result.root_input_type.columns,
            {
                'fastq1': wf_type.FILE,
                'fastq2': wf_type.FILE,
                'outbase': wf_type.STR,
                'ref': wf_type.FILE,
                'ref_fai': wf_type.FILE,
            },
        )
        self.assertTrue(isinstance(result.root_output_type, wf_type.RecordType))

    def test_batch_map_missing_field_reports_compile_time_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'batch_bad_field.swl', _BATCH_BAD_FIELD)
        result = Checker().load(path)
        self.assertTrue(any('Missing field on tab: nope' in err for err in result.errors))

    def test_batch_map_on_batch_workflow_reports_error(self):
        root = self._make_fixture_dir()
        self._write(root, 'batch_inner.swl', _BATCH_INNER)
        path = self._write(root, 'batch_on_batch.swl', _BATCH_ON_BATCH)
        result = Checker().load(path)
        self.assertTrue(any('map on batch workflow is not supported' in err for err in result.errors))

    def test_batch_map_requires_function_value(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'batch_non_function.swl', _BATCH_NON_FUNCTION)
        result = Checker().load(path)
        self.assertTrue(any('map requires a function value' in err for err in result.errors))

    def test_workflow_saturated_task_application_reports_not_a_function(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'apply_final.swl', _NON_FUNCTION_APPLY)
        result = Checker().load(path)
        self.assertIsNone(result.signature)
        self.assertIn('Workflow must evaluate to a function', result.errors)

    def test_simple_workflow_exposes_record_root_types(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'function.swl', _FUNCTION)
        result = Checker().load(path)
        self.assertTrue(isinstance(result.workflow_type, wf_type.FunctionType))
        self.assertTrue(isinstance(result.root_input_type, wf_type.RecordType))
        self.assertTrue(isinstance(result.root_output_type, wf_type.RecordType))

    def test_map_by_requires_existing_key_and_preserves_it(self):
        root = self._make_fixture_dir()
        self._write(root, 'group_ok.swl', 'align = import "align.sh"\n\\x ->\n    a = align x\n    { sample: x.sample, bam: a.bam }\n')
        path = self._write(root, 'map_by_ok.swl', 'g = import "group_ok.swl"\n\\xs ->\n    map_by g "sample" xs\n')
        result = Checker().load(path)
        self.assertEqual(result.errors, [])
        self.assertTrue(isinstance(result.root_input_type, wf_type.TableType))
        self.assertIn('sample', result.root_input_type.columns)

    def test_map_by_reports_missing_key_column(self):
        root = self._make_fixture_dir()
        self._write(root, 'group_nosample.sh', '# @ GroupNoSample\n# in\n#   bam [file]\n# out\n#   sample str = s\n#   bam file = out.bam\necho ok\n')
        self._write(root, 'group_nosample.swl', 'g = import "group_nosample.sh"\ng\n')
        path = self._write(root, 'map_by_missing_key.swl', 'g = import "group_nosample.swl"\n\\xs ->\n    map_by g "sample" xs\n')
        result = Checker().load(path)
        self.assertTrue(any('Missing field on tab: sample' in err for err in result.errors))

    def test_map_by_requires_preserving_grouping_key_in_output(self):
        root = self._make_fixture_dir()
        self._write(root, 'group_bad.sh', '# @ GroupBad\n# in\n#   sample [str]\n#   bam [file]\n# out\n#   bam file = out.bam\necho bad\n')
        self._write(root, 'group_bad.swl', 'g = import "group_bad.sh"\ng\n')
        path = self._write(root, 'map_by_missing_output_key.swl', 'g = import "group_bad.swl"\n\\xs ->\n    map_by g "sample" xs\n')
        result = Checker().load(path)
        self.assertTrue(any('map_by output must preserve grouping key: sample' in err for err in result.errors))

    # def test_table_update_reports_explicit_not_implemented(self):
    #     root = self._make_fixture_dir()
    #     path = self._write(root, 'tab_update.swl', 'align = import "align.sh"\n\\xs ->\n    ys = map align xs\n    ys // { extra: "x" }\n')
    #     result = Checker().load(path)
    #     self.assertTrue(any('table update semantics are not implemented' in err for err in result.errors))


if __name__ == '__main__':
    ut.main()
