from enum import Enum

TokenType = Enum('TokenType',
    [
        'id',
        'str',
        'num',
        'equal',
        'colon',
        'comma',
        'dot',
        'amp',
        'lparen',
        'rparen',
        'lbracket',
        'rbracket',
        'lbrace',
        'rbrace',
        'bslash',
        'arrow',
        'pipe',
        'indent',
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
        elif self.type == TokenType.indent:
            return f'indent({self.value})'
        elif self.type == TokenType.str:
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
        self.new_line = True

    def __iter__(self):
        return self

    def n_remaining(self):
        return len(self.s) - self.i

    def __next__(self):
        '''Obtain next token'''
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
                        return Token(TokenType.eol)

                # assess indentation level
                if self.new_line:
                    ws = self._until(lambda x: x == '\n' or not x.isspace())
                    if self.i < len(self.s) and self.s[self.i] == '\n':
                        # entire line is blank: advance and ignore line
                        self.i += 1
                        return self.__next__()
                    else:
                        # tab is equivalent to 4 spaces
                        ws = ws.replace('\t', '    ')
                        return Token(TokenType.indent, len(ws))

                # ignore whitespace in other contexts
                return self.__next__()

            self.new_line = False

            if s1 == '#':
                # ignore everything until end of line
                comment = self._jump(self.i, '\n')
                return self.__next__()

            if s1 == '&':
                return Token(TokenType.amp)
            
            if s1 == '(':
                return Token(TokenType.lparen)

            if s1 == ')':
                return Token(TokenType.rparen)
            
            if s1 == '{':
                return Token(TokenType.lbrace)

            if s1 == '}':
                return Token(TokenType.rbrace)
            
            if s1 == ':':
                return Token(TokenType.colon)

            if s1 == ',':
                return Token(TokenType.comma)

            if s1 == '.':
                return Token(TokenType.dot)

            if s1 == '\\':
                return Token(TokenType.bslash)

            if s1 == '=':
                return Token(TokenType.equal)

            # string literal
            if s1 == '"':
                value = self._match(self.i, '"')
                return Token(TokenType.str, value)

            # identifier
            if s1.isalpha():
                # advance index past next character that is not a valid
                # identifier character
                value = self._until(lambda x: not (x.isalnum() or x == '_'))
                return Token(TokenType.id, value)

            # tokens with 2 characters
            if self.n_remaining() >= 2:
                s2 = self.s[self.i:(self.i + 2)]
                # preemptively increment the index
                self.i += 1

                if s2 == '->':
                    return Token(TokenType.arrow)

                if s2 == '|>':
                    return Token(TokenType.pipe)
                
                # at this point, no 2-character match was found
                # move the index back to the previous state
                self.i -= 1

            raise Exception("Unrecognized character: '{}'".format(s1))

        elif self.i == len(self.s):
            self.i += 1
            return Token(TokenType.eof)

        else:
            raise StopIteration


    def _jump(self, start, r):
        '''Advance until substring r is found and return token not including r.
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

    def _match(self, start, r):
        '''Advance until substring r is found and return token not including r,
           returning exception if r is not found.
           Cursor position will be one past r.'''
        end = self.s.find(r, start)
        if end == -1:
            raise Exception('Unexpected EOF when searching for "{}"'.format(r))
        else:
            # advance index past the matching substring
            self.i = end + 1
            return self.s[start:end]

    def _until(self, predicate):
        '''Advance until predicate is true and return token,
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

