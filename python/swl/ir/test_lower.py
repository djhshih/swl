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

_COMPOSE = '''align = import "align.sh"
sub = import "sub.swl"
\\x ->
    a = align x
    b = sub (x // a)
    a // b
'''


class TestLower(ut.TestCase):
    def _files(self):
        root = os.path.abspath('/virtual')
        return {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sub.swl'): _SUB,
            os.path.join(root, 'main.swl'): _MAIN,
            os.path.join(root, 'chain.swl'): _CHAIN,
            os.path.join(root, 'compose.swl'): _COMPOSE,
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
        self.assertIs(body.bindings[0].value.function, body.result.fields['bam'].record.function)
        self.assertIs(body.bindings[1].value.function, body.result.fields['sbam'].record.function)

    def test_lower_chain_normalizes_to_compose(self):
        files, root = self._files()
        tree = parse_and_lower(_CHAIN, root, files)
        self.assertIsInstance(tree, ir.Compose)
        self.assertEqual([stage.name for stage in tree.stages], ['_s1', '_s2'])
        self.assertEqual(tree.param, '_input')

    def test_lower_explicit_compose_normalizes_to_compose(self):
        files, root = self._files()
        tree = parse_and_lower(_COMPOSE, root, files)
        self.assertIsInstance(tree, ir.Compose)
        self.assertEqual(tree.param, 'x')
        self.assertEqual([stage.name for stage in tree.stages], ['a', 'b'])


if __name__ == '__main__':
    ut.main()
