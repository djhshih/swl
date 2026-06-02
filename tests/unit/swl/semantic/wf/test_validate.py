import os
import tempfile
import unittest

from swl.semantic.wf.check import Checker
from swl.semantic.wf.validate import WorkflowInputValidationError, validate_workflow_inputs


_ALIGN = '''# @ Align
# in
#   fastq1, fastq2 file
#   ref file
#   ref_fai file
#   outbase str
# out
#   bam file = ${outbase}.bam
echo align
'''

_MERGE = '''# @ Merge
# in
#   bam [file]
#   outbase str
# out
#   bam file = ${outbase}.bam
echo merge
'''

_SIMPLE = '''align = import "align.sh"
\\x ->
    align x
'''

_BATCH = '''align = import "align.sh"
merge = import "merge.sh"
\\xs ->
    calls = map align xs
    merge { bam: calls.bam, outbase: "panel" }
'''


class TestWorkflowInputValidate(unittest.TestCase):
    def _fixture_dir(self):
        return tempfile.mkdtemp(prefix='swl-validate-')

    def _write(self, root, name, content):
        path = os.path.join(root, name)
        with open(path, 'w') as f:
            f.write(content)
        return path

    def test_validate_simple_workflow_inputs(self):
        root = self._fixture_dir()
        self._write(root, 'align.sh', _ALIGN)
        path = self._write(root, 'simple.swl', _SIMPLE)
        result = Checker().load(path)
        value = {
            'fastq1': 'a.fq',
            'fastq2': 'b.fq',
            'ref': 'hg38.fa',
            'ref_fai': 'hg38.fa.fai',
            'outbase': 'sample',
        }
        self.assertEqual(validate_workflow_inputs(result, value), value)

    def test_validate_batch_workflow_accepts_equal_length_columns(self):
        root = self._fixture_dir()
        self._write(root, 'align.sh', _ALIGN)
        self._write(root, 'merge.sh', _MERGE)
        path = self._write(root, 'batch.swl', _BATCH)
        result = Checker().load(path)
        value = {
            'fastq1': ['a_1.fq', 'b_1.fq'],
            'fastq2': ['a_2.fq', 'b_2.fq'],
            'ref': ['hg38.fa', 'hg38.fa'],
            'ref_fai': ['hg38.fa.fai', 'hg38.fa.fai'],
            'outbase': ['a', 'b'],
        }
        self.assertEqual(validate_workflow_inputs(result, value), value)

    def test_validate_batch_workflow_rejects_unequal_column_lengths(self):
        root = self._fixture_dir()
        self._write(root, 'align.sh', _ALIGN)
        self._write(root, 'merge.sh', _MERGE)
        path = self._write(root, 'batch.swl', _BATCH)
        result = Checker().load(path)
        with self.assertRaises(WorkflowInputValidationError) as ctx:
            validate_workflow_inputs(result, {
                'fastq1': ['a_1.fq', 'b_1.fq'],
                'fastq2': ['a_2.fq'],
                'ref': ['hg38.fa', 'hg38.fa'],
                'ref_fai': ['hg38.fa.fai', 'hg38.fa.fai'],
                'outbase': ['a', 'b'],
            })
        self.assertIn('equal length', str(ctx.exception))

    def test_validate_batch_workflow_rejects_non_array_column(self):
        root = self._fixture_dir()
        self._write(root, 'align.sh', _ALIGN)
        self._write(root, 'merge.sh', _MERGE)
        path = self._write(root, 'batch.swl', _BATCH)
        result = Checker().load(path)
        with self.assertRaises(WorkflowInputValidationError) as ctx:
            validate_workflow_inputs(result, {
                'fastq1': 'a_1.fq',
                'fastq2': ['a_2.fq'],
                'ref': ['hg38.fa'],
                'ref_fai': ['hg38.fa.fai'],
                'outbase': ['a'],
            })
        self.assertIn('array column', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
