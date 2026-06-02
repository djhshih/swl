import unittest as ut

from swl.syntax.wf.lexer import Lexer, Token, TokenType


class TestLexer(ut.TestCase):

    def test_assignment(self):
        lexer = Lexer('name1234 = "James"')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.id, 'name1234'),
            Token(TokenType.equal),
            Token(TokenType.str, 'James'),
            Token(TokenType.eof)]
        )

    def test_comment(self):
        lexer = Lexer('name  # ignored comment')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.id, 'name'),
            Token(TokenType.eof)]
        )

    def test_indented_comment_line(self):
        lexer = Lexer('\\x ->\n    # ignored\n    y')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.bslash),
            Token(TokenType.id, 'x'),
            Token(TokenType.arrow),
            Token(TokenType.bstart),
            Token(TokenType.id, 'y'),
            Token(TokenType.bend),
            Token(TokenType.eof)]
        )

    def test_hash_inside_string_is_not_comment(self):
        lexer = Lexer('x = "# not comment" # real comment')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.id, 'x'),
            Token(TokenType.equal),
            Token(TokenType.str, '# not comment'),
            Token(TokenType.eof)]
        )

    def test_escaped_quote_inside_string_is_not_terminator(self):
        lexer = Lexer('x = "a \\\"quoted\\\" # string" # real comment')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.id, 'x'),
            Token(TokenType.equal),
            Token(TokenType.str, 'a \\\"quoted\\\" # string'),
            Token(TokenType.eof)]
        )

    def test_function(self):
        lexer = Lexer('\\x ->\n    y = f x    x // y')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.bslash),
            Token(TokenType.id, 'x'),
            Token(TokenType.arrow),
            Token(TokenType.bstart),
            Token(TokenType.id, 'y'),
            Token(TokenType.equal),
            Token(TokenType.id, 'f'),
            Token(TokenType.id, 'x'),
            Token(TokenType.id, 'x'),
            Token(TokenType.update),
            Token(TokenType.id, 'y'),
            Token(TokenType.bend),
            Token(TokenType.eof)]
        )

    def test_chain(self):
        lexer = Lexer('a | b | c')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.id, 'a'),
            Token(TokenType.chain),
            Token(TokenType.id, 'b'),
            Token(TokenType.chain),
            Token(TokenType.id, 'c'),
            Token(TokenType.eof)]
        )

    def test_record(self):
        lexer = Lexer('{ align: { fastq1: "1.fq", fastq2: "2.fq" } }')
        self.assertEqual(
            [x for x in lexer],
            [Token(TokenType.lbrace),
            Token(TokenType.id, 'align'),
            Token(TokenType.colon),
            Token(TokenType.lbrace),
            Token(TokenType.id, 'fastq1'),
            Token(TokenType.colon),
            Token(TokenType.str, '1.fq'),
            Token(TokenType.comma),
            Token(TokenType.id, 'fastq2'),
            Token(TokenType.colon),
            Token(TokenType.str, '2.fq'),
            Token(TokenType.rbrace),
            Token(TokenType.rbrace),
            Token(TokenType.eof)]
        )


if __name__ == '__main__':
    ut.main()
