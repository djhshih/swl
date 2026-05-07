import os
import tempfile
import unittest as ut

from swl.semantic.wf.check import Checker


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
        return root

    def test_load_pipe_workflow(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'pipe.swl', _PIPE)
        result = Checker().load(path)
        self.assertEqual(sorted(result.imports.keys()), ['align', 'call', 'sort'])
        self.assertEqual(result.chain_errors, [])
        self.assertEqual(result.inferred_inputs, set())
        self.assertIn('bcf', result.signature.outputs)

    def test_load_function_workflow_infers_inputs(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'function.swl', _FUNCTION)
        result = Checker().load(path)
        self.assertEqual(result.chain_errors, [])
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
        self.assertEqual(result.chain_errors, [])
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
        self.assertEqual(result.chain_errors, [])
        self.assertIn('fastq1', result.signature.inputs)
        self.assertIn('fastq2', result.signature.inputs)
        self.assertIn('outbase', result.signature.inputs)
        self.assertNotIn('ref', result.signature.inputs)
        self.assertNotIn('ref_fai', result.signature.inputs)
        self.assertIn('bam', result.signature.outputs)

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
        self.assertEqual(result.chain_errors, [])
        self.assertIn('call', result.imports)
        self.assertEqual(result.imports['w'].kind, 'workflow')

    def test_direct_workflow_import_in_application(self):
        root = self._make_fixture_dir()
        self._write(root, 'task_result.swl', _TASK_RESULT_WORKFLOW)
        path = self._write(root, 'workflow_apply.swl', _WORKFLOW_APPLY)
        result = Checker().load(path)
        self.assertEqual(result.chain_errors, [])
        self.assertIn('fastq1', result.inferred_inputs)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)
        self.assertIn('bam', result.signature.outputs)

    def test_scalar_task_application_lifts_to_first_input(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'scalar_apply.swl', _SCALAR_APPLY)
        result = Checker().load(path)
        self.assertEqual(result.chain_errors, [])
        self.assertIsNotNone(result.signature)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)
        self.assertNotIn('fastq1', result.inferred_inputs)
        self.assertIn('bam', result.signature.outputs)

    def test_bad_pipe_reports_chain_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'bad_pipe.swl', _BAD_PIPE)
        result = Checker().load(path)
        self.assertEqual(len(result.chain_errors), 1)
        self.assertIn('outbase', result.chain_errors[0])

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

    def test_workflow_saturated_task_application_reports_not_a_function(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'apply_final.swl', _NON_FUNCTION_APPLY)
        result = Checker().load(path)
        self.assertIsNone(result.signature)
        self.assertIn('Workflow must evaluate to a function', result.errors)


if __name__ == '__main__':
    ut.main()
