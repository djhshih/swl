import unittest as ut

from swl.syntax.wf import node
from swl.syntax.wf.parser import Parser


class TestParser(ut.TestCase):

    def test_final_expr_in_block(self):
        src = "\\x ->\n    a = 1\n    { a: 1 }"
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_simple_chain(self):
        src = "align | sort"
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_simple_binding(self):
        src = "x = 1\ny"
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_final_binding_in_block_fails(self):
        src = "\\x ->\n    a = 1\n    b = 2"
        p = Parser()
        with self.assertRaises(ValueError):
            p.parse(src)

    def test_final_binding_in_outer_block_fails(self):
        src = "x = 1"
        p = Parser()
        with self.assertRaises(ValueError):
            p.parse(src)

    def test_leading_comment_parses(self):
        src = "# hello\nx = 1\ny"
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_trailing_comment_parses(self):
        src = "x = 1 # note\nx"
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_hash_inside_string_parses(self):
        src = 'x = "# not comment"\nx'
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_indented_comment_in_block_parses(self):
        src = "\\x ->\n    # note\n    y = x\n    y"
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_escaped_quote_and_hash_inside_string_parses(self):
        src = 'x = "a \\\"quoted\\\" # string"\nx'
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_trailing_comment_after_final_expr_parses(self):
        src = "x = 1\nx\n# done"
        p = Parser()
        result = p.parse(src)
        self.assertIsInstance(result, node.Expr)

    def test_final_binding_with_only_comment_after_still_fails(self):
        src = "x = 1\n# no final expr"
        p = Parser()
        with self.assertRaises(ValueError):
            p.parse(src)


if __name__ == '__main__':
    ut.main()
