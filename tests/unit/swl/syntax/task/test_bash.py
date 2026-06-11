import unittest as ut

from swl.syntax.task import bash, interpolation


class TestBashParser(ut.TestCase):

    def test_parse_assignments_and_commands(self):
        src = (
            'tmp=aligned.bam\n'
            'bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam\n'
        )
        script = bash.parse(src)
        self.assertEqual(len(script.statements), 2)
        self.assertEqual(script.statements[0].name, 'tmp')
        self.assertEqual(
            script.statements[0].value,
            interpolation.Word([interpolation.Literal('aligned.bam')])
        )
        self.assertEqual(
            script.statements[1].words,
            [
                interpolation.Word([interpolation.Var('cpu')]),
                interpolation.Word([interpolation.Var('ref')]),
                interpolation.Word([interpolation.Var('fastq1')]),
                interpolation.Word([interpolation.Var('fastq2')]),
                interpolation.Word([
                    interpolation.Var('outbase'),
                    interpolation.Literal('.bam')
                ]),
            ]
        )

    def test_parse_sort_body(self):
        src = (
            'samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam\n'
            'samtools index ${outbase}.bam ${outbase}.bai\n'
        )
        script = bash.parse(src)
        self.assertEqual(len(script.statements), 2)
        self.assertEqual(
            script.statements[0].words,
            [
                interpolation.Word([interpolation.Var('cpu')]),
                interpolation.Word([interpolation.Expr('memory / cpu')]),
                interpolation.Word([
                    interpolation.Var('outbase'),
                    interpolation.Literal('.bam')
                ]),
            ]
        )

    def test_empty_script(self):
        script = bash.parse('')
        self.assertEqual(len(script.statements), 0)

    def test_script_with_only_comment(self):
        script = bash.parse('# just a comment\n')
        self.assertEqual(len(script.statements), 0)


if __name__ == '__main__':
    ut.main()
