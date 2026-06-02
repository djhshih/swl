import unittest as ut

from swl.syntax.task.interpolation import Expr, Literal, Parser, Var, Word


class TestInterpolation(ut.TestCase):

    def test_var_with_suffix(self):
        result = Parser().parse_word('${outbase}.bam')
        self.assertEqual(result, Word([Var('outbase'), Literal('.bam')]))

    def test_simple_var(self):
        result = Parser().parse_word('${cpu}')
        self.assertEqual(result, Word([Var('cpu')]))

    def test_expr(self):
        result = Parser().parse_word('${memory / cpu}')
        self.assertEqual(result, Word([Expr('memory / cpu')]))

    def test_plain_literal(self):
        result = Parser().parse_word('aligned.bam')
        self.assertEqual(result, Word([Literal('aligned.bam')]))

    def test_unterminated_expr_fails(self):
        with self.assertRaises(ValueError):
            Parser().parse_word('${outbase')


if __name__ == '__main__':
    ut.main()
