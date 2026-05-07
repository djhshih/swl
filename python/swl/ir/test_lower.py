import os
import unittest as ut

from swl.ir import node as ir
from swl.ir.lower import parse_and_lower


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
    def _files(self):
        root = os.path.abspath('/virtual')
        return {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sub.swl'): _SUB,
            os.path.join(root, 'main.swl'): _MAIN,
            os.path.join(root, 'chain.swl'): _CHAIN,
        }, root

    def test_lower_lambda_workflow(self):
        files, root = self._files()
        tree = parse_and_lower(_MAIN, root, files)
        self.assertIsInstance(tree, ir.Lambda)
        self.assertEqual(tree.param, 'x')
        self.assertIsInstance(tree.body, ir.Block)

    def test_lower_imports_to_functions(self):
        files, root = self._files()
        tree = parse_and_lower(_MAIN, root, files)
        self.assertIsInstance(tree, ir.Lambda)
        body = tree.body
        self.assertIsInstance(body, ir.Block)
        self.assertEqual(len(body.bindings), 2)
        self.assertIsInstance(body.bindings[0].value, ir.Apply)
        self.assertIsInstance(body.bindings[0].value.function, ir.Function)
        self.assertEqual(body.bindings[0].value.function.kind, 'task')
        self.assertIsNone(body.bindings[0].value.function.body)
        self.assertIsInstance(body.bindings[1].value.function, ir.Function)
        self.assertEqual(body.bindings[1].value.function.kind, 'workflow')
        self.assertIsInstance(body.bindings[1].value.function.body, ir.Lambda)
        self.assertEqual(body.bindings[0].value.function.name, 'align')
        self.assertEqual(body.bindings[1].value.function.name, 'sub')

    def test_lower_chain_flattens(self):
        files, root = self._files()
        tree = parse_and_lower(_CHAIN, root, files)
        self.assertIsInstance(tree, ir.Chain)
        self.assertEqual(len(tree.items), 2)
        self.assertEqual(tree.items[0].kind, 'task')
        self.assertEqual(tree.items[1].kind, 'workflow')
        self.assertIsNone(tree.items[0].body)
        self.assertIsInstance(tree.items[1].body, ir.Lambda)


if __name__ == '__main__':
    ut.main()
