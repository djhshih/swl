from enum import Enum

class TokenType(Enum):
    eof = 1
    identifier = 2
    string = 3
    equal = 4
    colon = 5
    comma = 6
    lparen = 7
    rparen = 8
    lbrace = 9
    rbrace = 10
    function = 11
    arrow = 12
    comment = 14

class Token:
    def __init__(self, ttype, value=None):
        self.type = ttype
        self.value = value
    def __repr__(self):
        if self.type == TokenType.identifier:
            return f'id({self.value})'
        elif self.type == TokenType.string:
            return f'"{self.value}"'
        else:
            return f'{self.type.name}'

class Lexer:
    '''Iterator over tokens in a string'''

    def __init__(self, s, i = 0):
        '''Initialize lexer with string s'''
        self.s = s
        self.i = i

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

            # TODO parse whitespace for indent
            if s1.isspace():
                return self.__next__()

            if s1 == '#':
                # advance index past the end of line
                end = self.s.find('\n', self.i)
                if end == -1:
                    self.i = len(self.s)
                else:
                    self.i = end + 1
                return Token(TokenType.comment)
            
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

            if s1 == '\\':
                return Token(TokenType.function)

            if s1 == '=':
                return Token(TokenType.equal)

            # string literal
            if s1 == '"':
                # advance index past the matching quote
                end = self.s.find('"', self.i)
                if end == -1:
                    raise Exception('Unexpected EOF when searching for matching "')
                else:
                    start = self.i
                    self.i = end + 1
                    return Token(TokenType.string, self.s[start:end])

            # identifier
            if s1.isalpha():
                # advance index past next character that is not a valid
                # identifier character
                end = -1
                for j in range(self.i, len(self.s)):
                    sj = self.s[j]
                    if not (sj.isalnum() or sj == '_'):
                        end = j
                        break
                start = self.i - 1
                if end == -1:
                    self.i = len(self.s)
                    return Token(TokenType.identifier, self.s[start:])
                else:
                    self.i = end
                    return Token(TokenType.identifier, self.s[start:end])

            # tokens with 2 characters
            if self.n_remaining() >= 2:
                s2 = self.s[self.i:(self.i + 2)]
                # preemptively increment the index
                self.i += 1

                if s2 == '->':
                    return Token(TokenType.arrow)
                
                # at this point, no 2-character match was found
                # move the index back to the previous state
                self.i -= 1

            raise Exception("Unrecognized character: '{}'".format(s1))

        elif self.i == len(self.s):
            self.i += 1
            return Token(TokenType.eof)

        else:
            raise StopIteration

