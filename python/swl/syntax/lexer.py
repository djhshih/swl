from enum import Enum
import queue

TokenType = Enum('TokenType',
    [
        'id',
        'str',
        'num',
        'equal',
        'colon',
        'comma',
        'dot',
        'update',
        'lparen',
        'rparen',
        'lbracket',
        'rbracket',
        'lbrace',
        'rbrace',
        'bslash',
        'arrow',
        'pipe',
        'bstart',
        'bend',
        'eol',
        'eof'
    ]
)

# TODO add line and column number
class Token:
    def __init__(self, ttype, value=None):
        self.type = ttype
        self.value = value
    def __repr__(self):
        if self.type == TokenType.id:
            return f'id({self.value})'
        elif self.type == TokenType.str or self.type == TokenType.num:
            return f'"{self.value}"'
        else:
            return f'{self.type.name}'
    def __eq__(self, other):
        return self.type == other.type and self.value == other.value

class Lexer:
    '''Iterator over tokens in a string'''

    def __init__(self, s, i = 0):
        '''Initialize lexer with string s'''
        self.s = s
        self.i = i
        # indicate whether we are at a new line
        self.new_line = True
        # whether to ignore the next eol before a non-whitespace character
        self.ignore_eol = False
        # stack for keeping track of indentation level
        self.indent_stack = [0]
        # extra tokens
        self.tokens = queue.Queue()

    def __iter__(self):
        return self

    def size(self) -> int:
        '''Remaining size of input string'''
        return len(self.s) - self.i

    def __next__(self):
        '''Obtain next token'''

        # return queued tokens, if any
        if not self.tokens.empty():
            return self.tokens.get()

        if self.i < len(self.s):
            s1 = self.s[self.i]
            # self.i now refers to the next character
            self.i += 1

            if s1.isspace():

                # after new line character, set new line flag
                if s1 == '\n':
                    if self.new_line:
                        # ignore blank line
                        return self.__next__()
                    else:
                        self.new_line = True
                        if self.ignore_eol:
                            self.ignore_eol = False
                            return self.__next__()
                        else:
                            return Token(TokenType.eol)

                # assess indentation level
                if self.new_line:
                    ws = self._until(lambda x: x == '\n' or not x.isspace())
                    if self.i < len(self.s) and self.s[self.i] == '\n':
                        # entire line is blank: advance and ignore line
                        self.i += 1
                        return self.__next__()
                    else:
                        self.new_line = False
                        # tab is equivalent to 4 spaces
                        indent = len(ws.replace('\t', '    '))
                        if indent > self.indent_stack[-1]:
                            # add indent token since indentation level increased
                            self.indent_stack.append(indent)
                            return Token(TokenType.bstart)
                        elif indent < self.indent_stack[-1]:
                            # add dedent token since indentation level decreased
                            while indent < self.indent_stack[-1]:
                                self.indent_stack.pop()
                                self.tokens.put(Token(TokenType.bend))
                            return self.tokens.get()
                        else:
                            # indentation level remained the same
                            # advance to next token
                            return self.__next__()

                # ignore whitespace in other contexts
                return self.__next__()

            if s1 == '#':
                # ignore everything until end of line
                comment = self._jump(self.i, '\n')
                return self.__next__()

            # consider the situation where dedentation occured and the
            # first character of the new line is not a whitespace
            if self.new_line and self.indent_stack[-1] > 0:
                while 0 < self.indent_stack[-1]:
                    self.indent_stack.pop()
                    self.tokens.put(Token(TokenType.bend))
                self.i -= 1
                return self.tokens.get()

            # current character is not whitespace
            self.new_line = False
            self.ignore_eol = False

            # string literal
            if s1 == '"':
                value = self._match(self.i, '"')
                return Token(TokenType.str, value)

            # identifier
            if s1.isalpha() or s1 == '_':
                # advance index past next character that is not a valid
                # identifier character
                value = self._until(lambda x: not (x.isalnum() or x == '_'))
                return Token(TokenType.id, value)

            # number
            if s1.isdigit() or s1 == '-':
                value = self._until(lambda x: not x.isdigit())
                # '-' by itself is not a number
                if value != '-':
                    c = self.s[self.i]
                    fractional = ''
                    if c == '.':
                        # capture fractional part
                        self.i += 1
                        fractional = self._until(lambda x: not x.isdigit())
                        value += fractional
                    if c == 'e' or c == 'E':
                        # capture exponent part
                        self.i += 1
                        exponent = self._until(lambda x: not x.isdigit())
                        value += exponent
                    if fractional:
                        numeric = float(value)
                    else:
                        numeric = int(value)
                    return Token(TokenType.num, numeric)

            if s1 == '&':
                return Token(TokenType.update)
            
            if s1 == '(':
                return Token(TokenType.lparen)

            if s1 == ')':
                return Token(TokenType.rparen)
            
            if s1 == '{':
                self.ignore_eol = True
                return Token(TokenType.lbrace)

            if s1 == '}':
                return Token(TokenType.rbrace)
            
            if s1 == ':':
                return Token(TokenType.colon)

            if s1 == ',':
                self.ignore_eol = True
                return Token(TokenType.comma)

            if s1 == '.':
                return Token(TokenType.dot)

            if s1 == '\\':
                return Token(TokenType.bslash)

            if s1 == '=':
                return Token(TokenType.equal)

            if s1 == '|':
                self.ignore_eol = True
                return Token(TokenType.pipe)


            # tokens with 2 characters
            if self.size() >= 2:
                # self.i is already at the next position
                s2 = self.s[(self.i-1):(self.i + 1)]
                # preemptively increment the index
                self.i += 1

                if s2 == '->':
                    self.ignore_eol = True
                    return Token(TokenType.arrow)

                # at this point, no 2-character match was found
                # move the index back to the previous state
                self.i -= 1

            raise Exception("Unrecognized character: '{}'".format(s1))

        elif self.i == len(self.s):
            # we are at end of file, but we need to tear down before ending
            # add dedent tokens to bring indentation level to 0
            while 0 < self.indent_stack[-1]:
                self.indent_stack.pop()
                self.tokens.put(Token(TokenType.bend))
            # add end of file token
            self.tokens.put(Token(TokenType.eof))
            # increment so that we do not repeat the tear down
            self.i += 1
            return self.tokens.get()

        else:
            raise StopIteration


    def _jump(self, start: int, r: str) -> str:
        '''Advance until substring r is found and return substring not including r.
           Cursor position will be at r.'''
        end = self.s.find(r, start)
        if end == -1:
            # advance index past the end of file
            self.i = len(self.s)
            return self.s[start:]
        else:
            # advance index to the matching substring
            self.i = end
            return self.s[start:end]

    def _match(self, start: int, r: str) -> str:
        '''Advance until substring r is found and return substring not including r,
           returning exception if r is not found.
           Cursor position will be one past r.'''
        end = self.s.find(r, start)
        if end == -1:
            raise Exception('Unexpected EOF when searching for "{}"'.format(r))
        else:
            # advance index past the matching substring
            self.i = end + 1
            return self.s[start:end]

    def _until(self, predicate) -> str:
        '''Advance until predicate is true and return substring,
           assuming cursor i is at the next position.
           Cursor will be at the position where the predicate is true.'''
        end = -1
        for j in range(self.i, len(self.s)):
            sj = self.s[j]
            if predicate(sj):
                end = j
                break
        start = self.i - 1
        if end == -1:
            self.i = len(self.s)
            return self.s[start:]
        else:
            self.i = end
            return self.s[start:end]
        


import unittest as ut

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

    def test_function(self):
        lexer = Lexer('\\x ->\n    y = f x\n    x & y')
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
            Token(TokenType.eol),
            Token(TokenType.id, 'x'),
            Token(TokenType.update),
            Token(TokenType.id, 'y'),
            Token(TokenType.bend),
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

