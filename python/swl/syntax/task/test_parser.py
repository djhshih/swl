import os
import unittest as ut

from swl.syntax.task import node
from swl.syntax.task import interpolation
from swl.syntax.task.parser import Parser


_TESTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../../../tests')
)


class TestTaskParser(ut.TestCase):

    def test_parse_align_task(self):
        with open(os.path.join(_TESTS_DIR, 'align.sh')) as f:
            src = f.read()
        task = Parser().parse(src)

        self.assertEqual(task.annotation.doc, 'Align paired-end sequencing reads')
        self.assertEqual(len(task.annotation.sections), 3)
        self.assertTrue(task.body.startswith('bwa mem'))

        out_section = task.annotation.sections[1]
        self.assertEqual(out_section.kind, node.SectionType.OUT)
        self.assertEqual(len(out_section.params), 1)
        self.assertEqual(out_section.params[0].names, ['bam'])
        self.assertEqual(
            out_section.params[0].default,
            interpolation.Word([
                interpolation.Var('outbase'),
                interpolation.Literal('.bam'),
            ])
        )

    def test_multiline_description(self):
        src = (
            '#@ doc\n'
            '# in\n'
            '#   ref file\n'
            '#     | reference bwa index files\n'
            'echo hi\n'
        )
        task = Parser().parse(src)
        section = task.annotation.sections[0]
        self.assertEqual(section.params[0].desc, 'reference bwa index files')

    def test_missing_doc_fails(self):
        src = '# in\n#   bam file\ncmd\n'
        with self.assertRaises(ValueError):
            Parser().parse(src)


if __name__ == '__main__':
    ut.main()
