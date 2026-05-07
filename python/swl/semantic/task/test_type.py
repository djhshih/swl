import os
import unittest as ut

from swl.semantic.task.type import TypeChecker, TypeKind, signature_from_task
from swl.syntax.task.parser import Parser


_TESTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../../../tests')
)


class TestTaskType(ut.TestCase):

    def test_signature_from_align_task(self):
        with open(os.path.join(_TESTS_DIR, 'align.sh')) as f:
            task = Parser().parse(f.read())

        sig = signature_from_task(task)
        self.assertEqual(sig.inputs['fastq1'].type, TypeKind.FILE)
        self.assertEqual(sig.inputs['fastq2'].type, TypeKind.FILE)
        self.assertEqual(sig.inputs['outbase'].type, TypeKind.STR)
        self.assertEqual(sig.outputs['bam'].type, TypeKind.FILE)
        self.assertIn('cpu', sig.run)
        self.assertIsNotNone(sig.outputs['bam'].default)

    def test_duplicate_param_fails(self):
        src = (
            '#@ doc\n'
            '# in\n'
            '#   a file\n'
            '#   a file\n'
            'echo hi\n'
        )
        task = Parser().parse(src)
        with self.assertRaises(ValueError):
            signature_from_task(task)

    def test_type_checker_chain(self):
        parser = Parser()
        checker = TypeChecker()

        with open(os.path.join(_TESTS_DIR, 'align.sh')) as f:
            checker.add_task('align', signature_from_task(parser.parse(f.read())))
        with open(os.path.join(_TESTS_DIR, 'sort.sh')) as f:
            checker.add_task('sort', signature_from_task(parser.parse(f.read())))

        self.assertEqual(checker.check_chain('align', 'sort'), [])


if __name__ == '__main__':
    ut.main()
