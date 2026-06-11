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

    def test_literal_before_var(self):
        result = Parser().parse_word('prefix${var}')
        self.assertEqual(result, Word([Literal('prefix'), Var('var')]))

    def test_multiple_vars(self):
        result = Parser().parse_word('${a}${b}')
        self.assertEqual(result, Word([Var('a'), Var('b')]))

    def test_escaped_dollar_is_literal(self):
        result = Parser().parse_word('\\${notavar}')
        self.assertEqual(result, Word([Literal('\\'), Var('notavar')]))

    def test_empty_brace_fails(self):
        with self.assertRaises(ValueError):
            Parser().parse_word('${}')

    def test_complex_expression(self):
        result = Parser().parse_word('${a + b * c}')
        self.assertEqual(result, Word([Expr('a + b * c')]))

    def test_empty_string(self):
        result = Parser().parse_word('')
        self.assertEqual(result, Word([]))


if __name__ == '__main__':
    ut.main()
