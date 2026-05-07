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

_BAD_PIPE = '''align = import "bad_outbase.sh"
sort = import "sort.sh"
align | sort
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

    def test_bad_pipe_reports_chain_error(self):
        root = self._make_fixture_dir()
        path = self._write(root, 'bad_pipe.swl', _BAD_PIPE)
        result = Checker().load(path)
        self.assertEqual(len(result.chain_errors), 1)
        self.assertIn('outbase', result.chain_errors[0])


if __name__ == '__main__':
    ut.main()
