import os
import shutil
import tempfile
import unittest as ut

import swl.ir.node as ir
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


def _strip_ids(node, root_name=None):
    if isinstance(node, ir.Variable):
        return ('Variable', node.name, _strip_ids(node.value, root_name))
    if isinstance(node, ir.Ref):
        return ('Ref', node.name)
    if isinstance(node, ir.Lambda):
        return ('Lambda', '$root', _strip_ids(node.body, node.param))
    if isinstance(node, ir.Block):
        return ('Block', [_strip_ids(bind, root_name) for bind in node.bindings], _strip_ids(node.result, root_name))
    if isinstance(node, ir.Apply):
        return ('Apply', _strip_ids(node.function, root_name), _strip_ids(node.arg, root_name))
    if isinstance(node, ir.Update):
        return ('Update', _strip_ids(node.left, root_name), _strip_ids(node.right, root_name))
    if isinstance(node, ir.Record):
        return ('Record', {name: _strip_ids(value, root_name) for name, value in node.fields.items()})
    if isinstance(node, ir.Field):
        return ('Field', _strip_ids(node.record, root_name), node.name)
    if isinstance(node, ir.Function):
        return ('Function', node.name, node.kind)
    if isinstance(node, ir.Name):
        return ('Name', '$root' if node.name == root_name else node.name)
    if isinstance(node, ir.Input):
        return ('Input', node.name, node.type, node.desc)
    if isinstance(node, ir.Literal):
        return ('Literal', node.value)
    return node.__class__.__name__


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
        self.assertIsInstance(tree, ir.Block)
        inner = tree.result
        self.assertIsInstance(inner, ir.Lambda)
        self.assertEqual(inner.param, 'x')
        self.assertIsInstance(inner.body, ir.Block)

    def test_lower_imports_to_functions(self):
        files, root = self._files()
        tree = parse_and_lower(_MAIN, root, files)
        self.assertIsInstance(tree, ir.Block)
        self.assertEqual(len(tree.bindings), 2)
        self.assertIsInstance(tree.result, ir.Lambda)

    def test_lower_bindings_produce_variables_and_refs(self):
        files, root = self._files()
        tree = parse_and_lower(_MAIN, root, files)
        # Top-level import bindings wrap the main lambda
        self.assertIsInstance(tree.bindings[0], ir.Variable)
        self.assertIsInstance(tree.bindings[1], ir.Variable)
        inner = tree.result
        self.assertIsInstance(inner, ir.Lambda)
        body = inner.body
        self.assertIsInstance(body, ir.Block)
        self.assertIn('bam', body.result.fields)
        self.assertIn('sbam', body.result.fields)

    def test_lower_chain_normalizes_to_lambda_block(self):
        files, root = self._files()
        tree = parse_and_lower(_CHAIN, root, files)
        self.assertIsInstance(tree, ir.Block)
        inner = tree.result
        self.assertIsInstance(inner, ir.Lambda)
        self.assertEqual(inner.param, '_input')
        self.assertIsInstance(inner.body, ir.Block)
        self.assertEqual([bind.name for bind in inner.body.bindings], ['_s1', '_s2'])

    def test_lower_explicit_compose_normalizes_to_lambda_block(self):
        files, root = self._files()
        tree = parse_and_lower(_COMPOSE, root, files)
        self.assertIsInstance(tree, ir.Block)
        inner = tree.result
        self.assertIsInstance(inner, ir.Lambda)
        self.assertEqual(inner.param, 'x')
        self.assertIsInstance(inner.body, ir.Block)
        self.assertEqual([bind.name for bind in inner.body.bindings], ['a', 'b'])

    def test_chain_and_explicit_have_equivalent_lowered_shape(self):
        root = os.path.abspath('/virtual-equivalent')
        files = {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sort.sh'): '''# @ Sort
# in
#   bam file
#   outbase str
# out
#   bam file = ${outbase}.bam
#   bai file = ${outbase}.bai
echo sort
''',
            os.path.join(root, 'call.sh'): '''# @ Call
# in
#   bam file
#   ref file
#   ref_fai file
#   outbase str
# out
#   bcf file = ${outbase}.bcf
echo call
''',
        }
        chain = parse_and_lower('''align = import "align.sh"
sort = import "sort.sh"
call = import "call.sh"

align | sort | call
''', root, files)
        explicit = parse_and_lower('''align = import "align.sh"
sort = import "sort.sh"
call = import "call.sh"

\\x ->
    _s1 = align x
    _s2 = sort (x // _s1)
    _s3 = call (x // _s1 // _s2)
    _s1 // _s2 // _s3
''', root, files)
        self.assertEqual(_strip_ids(chain), _strip_ids(explicit))


    def test_gap1_import_shadowing_in_lambda_body(self):
        tmpdir = tempfile.mkdtemp()
        old_cwd = os.getcwd()
        try:
            os.chdir(tmpdir)
            task_b = '''# @ Sort
# in
#   bam file
# out
#   bam file = out.bam
echo sort
'''
            for name, content in [('a.sh', _ALIGN), ('b.sh', task_b)]:
                with open(os.path.join(tmpdir, name), 'w') as f:
                    f.write(content)
            src = 't = import "a.sh"\n\\x ->\n    t = import "b.sh"\n    t { bam: "x" }\n'
            tree = parse_and_lower(src, tmpdir, {})
            # Top-level: Block with outer t binding wrapping a Lambda
            self.assertIsInstance(tree, ir.Block)
            self.assertGreater(len(tree.bindings), 0)
            outer_binding = tree.bindings[0]
            self.assertIsInstance(outer_binding, ir.Variable)
            self.assertEqual(outer_binding.name, 't')
            self.assertIsInstance(outer_binding.value, ir.Function)
            self.assertEqual(outer_binding.value.name, 't')
            # Lambda body: inner t shadows with import "b.sh"
            block = tree.result
            self.assertIsInstance(block, ir.Lambda)
            body = block.body
            self.assertIsInstance(body, ir.Block)
            inner_binding = body.bindings[0]
            self.assertIsInstance(inner_binding, ir.Variable)
            self.assertEqual(inner_binding.name, 't')
            self.assertIsInstance(inner_binding.value, ir.Function)
            self.assertEqual(inner_binding.value.name, 'b')
        finally:
            os.chdir(old_cwd)
            shutil.rmtree(tmpdir)


if __name__ == '__main__':
    ut.main()
