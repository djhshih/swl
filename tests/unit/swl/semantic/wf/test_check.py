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


    # ====================================================================
    # Category 1: Basic binding and scope isolation
    # ====================================================================

    def test_1_1_single_binding_resolves(self):
        r = Checker().load_content('x = 1\n\\y -> y', '/tmp/test_1_1.swl')
        self.assertEqual(r.errors, [])

    def test_1_2_multiple_bindings_same_scope(self):
        r = Checker().load_content('x = 1\ny = 2\n\\z -> z', '/tmp/test_1_2.swl')
        self.assertEqual(r.errors, [])

    def test_1_3_binding_refers_to_prior_binding(self):
        r = Checker().load_content('x = 1\ny = x\n\\z -> z', '/tmp/test_1_3.swl')
        self.assertEqual(r.errors, [])

    def test_1_4_binding_refers_to_outer_scope(self):
        r = Checker().load_content('x = 1\n\\y -> x', '/tmp/test_1_4.swl')
        self.assertEqual(r.errors, [])

    def test_1_5_final_expr_is_binding_raises_parse_error(self):
        with self.assertRaises(ValueError) as ctx:
            Checker().load_content('x = 1\ny = 2', '/tmp/test_1_5.swl')
        self.assertIn('Final line in a block must be an expr', str(ctx.exception))

    # ====================================================================
    # Category 2: Duplicate detection
    # ====================================================================

    def test_2_1_duplicate_binding(self):
        r = Checker().load_content('x = 1\nx = 2\n\\y -> y', '/tmp/test_2_1.swl')
        self.assertIn('Duplicate binding in scope: x', r.errors)

    def test_2_3_duplicate_binding_import_then_scalar(self):
        root = tempfile.mkdtemp()
        self._write(root, 'a.sh', _ALIGN)
        path = self._write(root, 'test_2_3.swl', 'x = import "a.sh"\nx = 1\n\\y -> y')
        r = Checker().load(path)
        self.assertIn('Duplicate binding in scope: x', r.errors)

    def test_2_4_block_inside_lambda_shadowing_allowed(self):
        r = Checker().load_content('\\x ->\n    x = 1\n    x', '/tmp/test_2_4.swl')
        self.assertEqual(r.errors, [])

    # ====================================================================
    # Category 3: Shadowing
    # ====================================================================

    def test_3_1_lambda_param_shadows_outer_binding(self):
        r = Checker().load_content('x = 1\n\\x -> x', '/tmp/test_3_1.swl')
        self.assertEqual(r.errors, [])

    def test_3_2_lambda_body_binding_shadows_outer_binding(self):
        r = Checker().load_content('x = 1\n\\y ->\n    x = 2\n    x', '/tmp/test_3_2.swl')
        self.assertEqual(r.errors, [])

    def test_3_3_deeply_nested_shadowing(self):
        r = Checker().load_content('a = 1\n\\b ->\n    a = 2\n    \\c ->\n        a = 3\n        a', '/tmp/test_3_3.swl')
        self.assertEqual(r.errors, [])

    def test_3_4_chain_desugaring_does_not_leak_generated_names(self):
        root = tempfile.mkdtemp()
        self._write(root, 'chainable.sh',
                    '# @ Chainable\n# in\n#   x file\n# out\n#   y file = out.y\necho chain')
        path = self._write(root, 'test_3_4.swl',
                           'a = import "chainable.sh"\nb = import "chainable.sh"\na | b')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_3_5_block_inside_lambda_shadows_outer(self):
        r = Checker().load_content('\\x ->\n    x = 42\n    x', '/tmp/test_3_5.swl')
        self.assertEqual(r.errors, [])

    # ====================================================================
    # Category 4: Import scope
    # ====================================================================

    def _scope_fixture_dir(self):
        td = tempfile.TemporaryDirectory()
        self.addCleanup(td.cleanup)
        root = td.name
        self._write(root, 'a.sh',
                    '# @ A\n# in\n#   x file\n# out\n#   y file = out.y\necho a')
        self._write(root, 'b.sh',
                    '# @ B\n# in\n#   x file\n# out\n#   y file = out.y\necho b')
        return root

    def test_4_1_import_inside_lambda_new_name(self):
        root = self._scope_fixture_dir()
        path = self._write(root, 'test_4_1.swl', '\\x ->\n    t = import "a.sh"\n    t x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_4_2_import_inside_lambda_direct_apply(self):
        root = self._scope_fixture_dir()
        path = self._write(root, 'test_4_2.swl', '\\x -> import "a.sh" x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_4_3_import_inside_lambda_shadows_top_level_import(self):
        root = self._scope_fixture_dir()
        path = self._write(root, 'test_4_3.swl',
                           't = import "a.sh"\n\\x ->\n    t = import "b.sh"\n    t x')
        r = Checker().load(path)
        # Scope allows shadowing; lowerer has a bug (Gap 1) but checker should pass
        self.assertNotIn('Duplicate binding in scope: t', r.errors)

    def test_4_4_import_inside_lambda_shadows_nonimport_outer(self):
        root = self._scope_fixture_dir()
        path = self._write(root, 'test_4_4.swl',
                           'x = 42\n\\y ->\n    x = import "a.sh"\n    x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_4_5_multiple_imports_same_file_different_names(self):
        root = self._scope_fixture_dir()
        path = self._write(root, 'test_4_5.swl',
                           'a = import "a.sh"\nb = import "a.sh"\n\\x -> x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_4_6_import_inside_lambda_same_name_as_top_level_nonimport(self):
        root = self._scope_fixture_dir()
        path = self._write(root, 'test_4_6.swl',
                           'x = 42\n\\y ->\n    x = import "a.sh"\n    x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    # ====================================================================
    # Category 5: Forward and self references
    # ====================================================================

    def test_5_1_forward_reference_across_bindings(self):
        r = Checker().load_content('x = y\ny = 1\n\\z -> z', '/tmp/test_5_1.swl')
        self.assertTrue(
            any('Forward reference' in err for err in r.errors),
        )

    def test_5_2_self_reference(self):
        r = Checker().load_content('x = x\n\\y -> y', '/tmp/test_5_2.swl')
        self.assertTrue(
            any('references itself' in err for err in r.errors),
        )

    def test_5_3_valid_cross_reference(self):
        r = Checker().load_content('x = 1\ny = x\n\\z -> z', '/tmp/test_5_3.swl')
        self.assertEqual(r.errors, [])

    def test_5_4_forward_reference_inside_lambda_body(self):
        r = Checker().load_content('\\a ->\n    x = y\n    y = 1\n    x', '/tmp/test_5_4.swl')
        self.assertTrue(
            any('Forward reference' in err for err in r.errors),
        )

    def test_5_5_self_reference_inside_lambda_body(self):
        r = Checker().load_content('\\a ->\n    x = x\n    x', '/tmp/test_5_5.swl')
        self.assertTrue(
            any('references itself' in err for err in r.errors),
        )

    # ====================================================================
    # Category 6: Annotation language scope (flat single scope, no shadowing)
    # ====================================================================

    def test_6_1_duplicate_name_in_in_section(self):
        root = tempfile.mkdtemp()
        self._write(root, 'bad.sh',
                    '# @ Bad\n# in\n#   x file\n#   x str\necho test')
        path = self._write(root, 'test_6_1.swl', 't = import "bad.sh"\n\\x -> x')
        with self.assertRaises(ValueError) as ctx:
            Checker().load(path)
        self.assertIn('Duplicate input parameter: x', str(ctx.exception))

    def test_6_2_cross_section_in_out_allowed(self):
        root = tempfile.mkdtemp()
        self._write(root, 'bad.sh',
                    '# @ Bad\n# in\n#   x file\n# out\n#   x file = out.txt\necho test')
        path = self._write(root, 'test_6_2.swl', 't = import "bad.sh"\n\\x -> x')
        # Allowed: many tasks share names between input and output
        # (e.g. sort takes bam, produces bam). Only in/run and out/run
        # cross-section duplicates are errors.
        Checker().load(path)

    def test_6_3_cross_section_in_run_duplicate_must_error(self):
        root = tempfile.mkdtemp()
        self._write(root, 'bad.sh',
                    '# @ Bad\n# in\n#   x file\n# run\n#   x = 2\necho test')
        path = self._write(root, 'test_6_3.swl', 't = import "bad.sh"\n\\x -> x')
        with self.assertRaises(ValueError):
            Checker().load(path)

    def test_6_4_cross_section_out_run_duplicate_must_error(self):
        root = tempfile.mkdtemp()
        self._write(root, 'bad.sh',
                    '# @ Bad\n# out\n#   x file = out.txt\n# run\n#   x = 2\necho test')
        path = self._write(root, 'test_6_4.swl', 't = import "bad.sh"\n\\x -> x')
        with self.assertRaises(ValueError):
            Checker().load(path)

    def test_6_5_duplicate_in_run_section_errors(self):
        root = tempfile.mkdtemp()
        self._write(root, 'bad.sh',
                    '# @ Bad\n# run\n#   cpu = 2\n#   cpu = 4\necho test\necho $cpu')
        path = self._write(root, 'test_6_5.swl', 't = import "bad.sh"\n\\x -> x')
        with self.assertRaises(ValueError) as ctx:
            Checker().load(path)
        self.assertIn('Duplicate run parameter: cpu', str(ctx.exception))

    def test_6_6_missing_type_on_in_param(self):
        root = tempfile.mkdtemp()
        self._write(root, 'bad.sh',
                    '# @ Bad\n# in\n#   x\necho test')
        path = self._write(root, 'test_6_6.swl', 't = import "bad.sh"\n\\x -> x')
        with self.assertRaises(ValueError) as ctx:
            Checker().load(path)
        self.assertIn('must have a type annotation', str(ctx.exception))

    def test_6_7_run_name_duplicates_in_must_error(self):
        root = tempfile.mkdtemp()
        self._write(root, 'bad.sh',
                    '# @ Bad\n# in\n#   cpu int\n# run\n#   cpu = 4\necho test')
        path = self._write(root, 'test_6_7.swl', 't = import "bad.sh"\n\\x -> x')
        with self.assertRaises(ValueError):
            Checker().load(path)

    # ====================================================================
    # Category 7: Interpolation scope
    # ====================================================================

    def test_7_1_output_default_references_input_param(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# in\n#   outbase str\n# out\n#   bam file = ${outbase}.bam\necho test')
        path = self._write(root, 'test_7_1.swl', 't = import "task.sh"\n\\x -> x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_7_2_output_default_references_run_param(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# run\n#   cpu = 2\n# out\n#   log file = out_${cpu}.log\necho test')
        path = self._write(root, 'test_7_2.swl', 't = import "task.sh"\n\\x -> x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_7_3_output_default_references_nonexistent_var(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# out\n#   bam file = ${nonexistent}.bam\necho test')
        path = self._write(root, 'test_7_3.swl', 't = import "task.sh"\n\\x -> x')
        with self.assertRaises(ValueError) as ctx:
            Checker().load(path)
        self.assertIn('Unresolved variable', str(ctx.exception))

    def test_7_4_command_block_references_assignment_lhs(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# out\n#   result file = out.txt\nx=foo\necho ${x}')
        path = self._write(root, 'test_7_4.swl', 't = import "task.sh"\n\\x -> x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_7_5_command_block_references_shell_builtin(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# out\n#   result file = out.txt\necho ${HOME}')
        path = self._write(root, 'test_7_5.swl', 't = import "task.sh"\n\\x -> x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_7_6_command_block_references_unknown_var(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# out\n#   result file = out.txt\necho ${nonexistent}')
        path = self._write(root, 'test_7_6.swl', 't = import "task.sh"\n\\x -> x')
        with self.assertRaises(ValueError) as ctx:
            Checker().load(path)
        self.assertIn('Unresolved variable', str(ctx.exception))

    def test_7_7_run_param_default_references_input(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# in\n#   threads int\n# run\n#   cpu = ${threads}\necho test')
        path = self._write(root, 'test_7_7.swl', 't = import "task.sh"\n\\x -> x')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    # ====================================================================
    # Category 8: Workflow well-formedness (final expr must be a function)
    # ====================================================================
    # Gap 2: Run-param and output default variable validation
    # ====================================================================

    def test_run_param_default_nonexistent_var(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# run\n#   cpu = ${nonexistent}\necho test')
        path = self._write(root, 'test_run_param_default_nonexistent_var.swl',
                           't = import "task.sh"\n\\x -> x')
        with self.assertRaises(ValueError) as ctx:
            Checker().load(path)
        self.assertIn('Unresolved variable', str(ctx.exception))

    def test_output_default_nonexistent_var_errors(self):
        root = tempfile.mkdtemp()
        self._write(root, 'task.sh',
                    '# @ Test\n# out\n#   bam file = ${nonexistent}.bam\necho test')
        path = self._write(root, 'test_output_default_nonexistent_var_errors.swl',
                           't = import "task.sh"\n\\x -> x')
        with self.assertRaises(ValueError) as ctx:
            Checker().load(path)
        self.assertIn('Unresolved variable', str(ctx.exception))

    # ====================================================================

    def test_8_1_scalar_as_final_expression(self):
        r = Checker().load_content('42', '/tmp/test_8_1.swl')
        self.assertIn('Workflow must evaluate to a function', r.errors)

    def test_8_2_record_as_final_expression(self):
        r = Checker().load_content('{ foo: 1 }', '/tmp/test_8_2.swl')
        self.assertIn('Workflow must evaluate to a function', r.errors)

    def test_8_3_saturated_task_application_as_final_expression(self):
        root = tempfile.mkdtemp()
        self._write(root, 'align.sh', _ALIGN)
        path = self._write(root, 'test_8_3.swl',
                           't = import "align.sh"\nt { fastq1: "a", fastq2: "b", ref: "r", ref_fai: "r.fai", outbase: "o" }')
        r = Checker().load(path)
        self.assertIsNone(r.signature)
        self.assertIn('Workflow must evaluate to a function', r.errors)

    def test_8_4_lambda_as_final_expression(self):
        r = Checker().load_content('\\x -> x', '/tmp/test_8_4.swl')
        self.assertIsNotNone(r.signature)
        self.assertEqual(r.errors, [])

    def test_8_5_imported_workflow_as_final_expression(self):
        root = tempfile.mkdtemp()
        self._write(root, 'sub.swl', '\\x ->\n    { y: x.a }')
        path = self._write(root, 'test_8_5.swl', 'w = import "sub.swl"\nw')
        r = Checker().load(path)
        self.assertIsNotNone(r.signature)
        self.assertEqual(r.errors, [])

    def test_8_6_partially_applied_task_as_final_expression(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'test_8_6.swl',
                           't = import "align.sh"\nt { ref: "r", ref_fai: "r.fai" }')
        r = Checker().load(path)
        self.assertIsNotNone(r.signature)
        self.assertEqual(r.errors, [])

    # ====================================================================
    # Category 9: Edge cases with chains and updates
    # ====================================================================

    def test_9_1_chain_with_non_imported_tasks_evaluates(self):
        root = tempfile.mkdtemp()
        self._write(root, 'a.sh',
                    '# @ A\n# in\n#   x file\n# out\n#   y file = out.y\necho a')
        self._write(root, 'b.sh',
                    '# @ B\n# in\n#   x file\n# out\n#   y file = out.y\necho b')
        path = self._write(root, 'test_9_1.swl',
                           'a = import "a.sh"\nb = import "b.sh"\na | b')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_9_2_chain_respects_binding_order(self):
        root = tempfile.mkdtemp()
        self._write(root, 'a.sh',
                    '# @ A\n# in\n#   x file\n# out\n#   y file = out.y\necho a')
        self._write(root, 'b.sh',
                    '# @ B\n# in\n#   x file\n# out\n#   y file = out.y\necho b')
        self._write(root, 'c.sh',
                    '# @ C\n# in\n#   x file\n# out\n#   y file = out.y\necho c')
        path = self._write(root, 'test_9_2.swl',
                           'a = import "a.sh"\nb = import "b.sh"\nc = import "c.sh"\na | b | c')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    # ====================================================================
    # Category 10: Type-level scope
    # ====================================================================

    def test_10_1_chain_type_check_compatible(self):
        root = tempfile.mkdtemp()
        self._write(root, 'a.sh',
                    '# @ A\n# in\n#   x file\n# out\n#   y file = out.y\necho a')
        self._write(root, 'b.sh',
                    '# @ B\n# in\n#   x file\n# out\n#   y file = out.y\necho b')
        path = self._write(root, 'test_10_1.swl',
                           'a = import "a.sh"\nb = import "b.sh"\na | b')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    def test_10_2_chain_type_check_mismatch(self):
        root = tempfile.mkdtemp()
        self._write(root, 'a.sh',
                    '# @ A\n# in\n#   x file\n# out\n#   x int = 42\necho a')
        self._write(root, 'b.sh',
                    '# @ B\n# in\n#   x file\n# out\n#   y file = out.y\necho b')
        path = self._write(root, 'test_10_2.swl',
                           'a = import "a.sh"\nb = import "b.sh"\na | b')
        r = Checker().load(path)
        self.assertTrue(
            any('type' in err.lower() for err in r.errors),
            f'Expected type error, got: {r.errors}',
        )

    def test_10_3_type_resolution_for_imported_task_chain(self):
        root = tempfile.mkdtemp()
        self._write(root, 'align.sh', _ALIGN)
        self._write(root, 'sort.sh', _SORT)
        path = self._write(root, 'test_10_3.swl',
                           'align = import "align.sh"\nsort = import "sort.sh"\nalign | sort')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])

    # ====================================================================
    # Gap 3: import as standalone expression (not bound to a name)
    # ====================================================================

    def test_gap3_standalone_import_as_final_expression(self):
        root = tempfile.mkdtemp()
        self._write(root, 'align.sh', _ALIGN)
        path = self._write(root, 'test_gap3.swl',
                           'import "align.sh"')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])
        self.assertIsNotNone(r.signature)

    def test_gap3_standalone_import_applied_directly(self):
        root = tempfile.mkdtemp()
        self._write(root, 'align.sh', _ALIGN)
        path = self._write(root, 'test_gap3b.swl',
                           '\\x -> import "align.sh" { fastq1: x }')
        r = Checker().load(path)
        self.assertEqual(r.errors, [])


if __name__ == '__main__':
    ut.main()
