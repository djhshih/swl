import os
import tempfile
import unittest as ut

from swl.ir import node as ir
from swl.ir.lower import lower_file


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

_SUB = '''align = import "align.sh"
\\x ->
    a = align x
    { bam: a.bam }
'''

_MAIN = '''align = import "align.sh"
sub = import "sub.swl"
\\x ->
    a = align x
    b = sub x
    { bam: a.bam, sbam: b.bam }
'''

_CHAIN = '''align = import "align.sh"
sub = import "sub.swl"
align | sub
'''


class TestLower(ut.TestCase):
    def _write(self, root, name, content):
        path = os.path.join(root, name)
        with open(path, 'w') as f:
            f.write(content)
        return path

    def _fixture(self):
        td = tempfile.TemporaryDirectory()
        self.addCleanup(td.cleanup)
        root = td.name
        self._write(root, 'align.sh', _ALIGN)
        self._write(root, 'sub.swl', _SUB)
        return root

    def test_lower_lambda_workflow(self):
        root = self._fixture()
        path = self._write(root, 'main.swl', _MAIN)
        tree = lower_file(path)
        self.assertIsInstance(tree, ir.Block)
        self.assertEqual(len(tree.bindings), 2)
        self.assertIsInstance(tree.result, ir.Lambda)
        self.assertIsInstance(tree.result.body, ir.Block)

    def test_lower_imports_to_single_import_node(self):
        root = self._fixture()
        path = self._write(root, 'main.swl', _MAIN)
        tree = lower_file(path)
        self.assertIsInstance(tree.bindings[0].value, ir.Import)
        self.assertIsInstance(tree.bindings[1].value, ir.Import)
        body = tree.result.body
        self.assertIsInstance(body.bindings[0].value, ir.Apply)
        self.assertIsInstance(body.bindings[0].value.function, ir.Import)
        self.assertEqual(body.bindings[0].value.function.kind, 'task')
        self.assertIsInstance(body.bindings[1].value.function, ir.Import)
        self.assertEqual(body.bindings[1].value.function.kind, 'workflow')
        self.assertEqual(body.bindings[0].value.function.name, 'align')
        self.assertEqual(body.bindings[1].value.function.name, 'sub')

    def test_lower_chain_flattens(self):
        root = self._fixture()
        path = self._write(root, 'chain.swl', _CHAIN)
        tree = lower_file(path)
        self.assertIsInstance(tree.result, ir.Chain)
        self.assertEqual(len(tree.result.items), 2)
        self.assertEqual(tree.result.items[0].kind, 'task')
        self.assertEqual(tree.result.items[1].kind, 'workflow')


if __name__ == '__main__':
    ut.main()
