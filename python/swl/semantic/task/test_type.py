import unittest as ut

from swl.semantic.task.type import TypeChecker, TypeKind, signature_from_task
from swl.syntax.task.parser import Parser


_ALIGN = (
    '# @ Align\n'
    '# in\n'
    '#   fastq1, fastq2 file\n'
    '#   outbase str\n'
    '# out\n'
    '#   bam file = ${outbase}.bam\n'
    '# run\n'
    '#   cpu = 2\n'
    'echo hi\n'
)

_SORT = (
    '# @ Sort\n'
    '# in\n'
    '#   bam file\n'
    '#   outbase str\n'
    '# out\n'
    '#   bam file = ${outbase}.bam\n'
    'echo hi\n'
)


class TestTaskType(ut.TestCase):

    def test_signature_from_align_task(self):
        task = Parser().parse(_ALIGN)
        sig = signature_from_task(task)
        self.assertEqual(sig.inputs['fastq1'].type, TypeKind.FILE)
        self.assertEqual(sig.inputs['fastq2'].type, TypeKind.FILE)
        self.assertEqual(sig.inputs['outbase'].type, TypeKind.STR)
        self.assertEqual(sig.outputs['bam'].type, TypeKind.FILE)
        self.assertIn('cpu', sig.run)
        self.assertIsNotNone(sig.outputs['bam'].default)

    def test_duplicate_input_param_fails(self):
        src = (
            '#@ doc\n'
            '# in\n'
            '#   a file\n'
            '#   a file\n'
            'echo hi\n'
        )
        task = Parser().parse(src)
        with self.assertRaises(ValueError) as cm:
            signature_from_task(task)
        self.assertEqual(str(cm.exception), 'Duplicate input parameter: a')

    def test_duplicate_output_param_fails(self):
        src = (
            '#@ doc\n'
            '# out\n'
            '#   bam file = result.bam\n'
            '#   bam file = result2.bam\n'
            'echo hi\n'
        )
        task = Parser().parse(src)
        with self.assertRaises(ValueError) as cm:
            signature_from_task(task)
        self.assertEqual(str(cm.exception), 'Duplicate output parameter: bam')

    def test_input_and_output_may_share_name(self):
        src = (
            '#@ doc\n'
            '# in\n'
            '#   bam file\n'
            '# out\n'
            '#   bam file = result.bam\n'
            'echo hi\n'
        )
        task = Parser().parse(src)
        sig = signature_from_task(task)
        self.assertIn('bam', sig.inputs)
        self.assertIn('bam', sig.outputs)

    def test_output_without_default_fails(self):
        src = (
            '#@ doc\n'
            '# out\n'
            '#   bam file\n'
            'echo hi\n'
        )
        task = Parser().parse(src)
        with self.assertRaises(ValueError) as cm:
            signature_from_task(task)
        self.assertEqual(str(cm.exception), 'Output parameter must have a default: bam')

    def test_type_checker_chain(self):
        parser = Parser()
        checker = TypeChecker()
        checker.add_task('align', signature_from_task(parser.parse(_ALIGN)))
        checker.add_task('sort', signature_from_task(parser.parse(_SORT)))
        self.assertEqual(checker.check_chain('align', 'sort'), [])


if __name__ == '__main__':
    ut.main()
