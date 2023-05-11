from swl.parse import node
from swl.parse import lexer as lex
from swl.parse.lexer import Token, TokenType

class Parser:

    def parse(self, s: str) -> node.Expr:
        '''Parse string `s` into an abstract syntax tree.'''
        self.lexer = lex.Lexer(s)
        self.queue = [next(self.lexer)]
        return self._parse_block()

    def eof(self) -> bool:
        '''Return whether parser is at end of file.'''
        if not self.queue: return True
        return self.queue[0].type == TokenType.eof

    def _at(self, i = 0) -> Token:
        '''Peek token at position `i`.'''
        if i < len(self.queue):
            # return extracted token
            return self.queue[i]
        else:
            # extract tokens until target position
            while len(self.queue) <= i:
                token = next(self.lexer)
                self.queue.append(token) 
            return token
    
    def _eat(self) -> Token:
        '''Advance and return next token.'''
        if self.queue:
            # eat previously extracted token if available
            token = self.queue.pop(0)
        else:
            # otherwise, eat next token
            token = next(self.lexer)

        # ensure that queue is not empty unless we are at eof
        if not self.queue:
            try:
                self.queue.append(next(self.lexer))
            except StopIteration:
                pass

        return token

    def _optional(self, token_type: TokenType) -> Token:
        '''Eat token if it has the expected type'''
        if self._at().type == token_type:
            return self._eat()
        return None

    def _expect(self, token_type: TokenType) -> Token:
        '''Eat a token and check that it has the expected type'''
        token = self._eat()
        if not token:
            raise ValueError(f'Parsing failed due to undefined token')
        if token.type != token_type:
            raise ValueError(
                f'Parsing failed due to unexpected token: expected token type {token_type.name}, but got token {token}'
            )
        return token

    def _until(self, token_type: TokenType) -> bool:
        '''Return whether current token is of type `token_type`.'''
        return not self.eof() and self.queue[0].type != token_type

    def _find(self, token_type: TokenType) -> int:
        '''Find and return index of token with type `token_type`,
           but do not advance.'''
        i = 0
        while not self.eof():
            token = self._at(i)
            if token.type == token_type:
                return i
            i += 1
        return -1

    def _parse_block(self, inner=False):
        '''Parse a block of code, which can be outer-level or inner-level.'''
        exprs = []
        if inner:
            # check for block start token
            self._expect(TokenType.bstart)
        while True:
            exprs.append(self._parse_expr())
            if self.eof(): break
            # eol right before eof is optional
            self._expect(TokenType.eol)
        if inner:
            # check for block end token
            self._exepct(TokenType.bend)
        return node.Block(exprs)

    def _parse_expr(self):
        t = self.queue[0].type

        if t == TokenType.bslash:
            return self._parse_function()

        elif t == TokenType.id and self._at(1).type == TokenType.equal:
            return self._parse_binding()

        else:
            return self._parse_simple_expr()

    # order of precedence for operators
    # from outer-most parser (lowest precedence)
    # to inner-most parser (highest precedence)
    # pipe, update, apply, get

    def _parse_id(self):
        iden = self._expect(TokenType.id)
        return node.Identifier(iden.value)

    def _parse_record(self):
        self._eat()  # eat open brace
        d = {}
        while self._until(TokenType.rbrace):
            iden = self._parse_id()
            self._expect(TokenType.colon)
            value = self._parse_expr()
            # comma is required after each key-value pair
            # unless we just parsed the final key-value pair
            if self._at().type != TokenType.rbrace:
                self._expect(TokenType.comma)
            d[iden.name] = value
        self._expect(TokenType.rbrace)
        return node.Record(d)

    def _parse_group(self):
        self._eat()  # eat open parenthesis
        value = self._parse_expr()
        self._expect(TokenType.rparen)
        return value

    def _parse_function(self):
        self._eat()  # eat backslash
        param = self._parse_id()
        self._expect(TokenType.arrow)
        block = self._parse_block(True)
        return node.Function(iden, block)

    def _parse_binding(self):
        iden = self._parse_id()
        self._expect(TokenType.equal)
        expr = self._parse_expr()
        return node.Binding(iden, expr)

    def _parse_pipe_expr(self):
        # TODO
        return self._parse_update_expr()

    def _parse_update_expr(self):
        # TODO
        return self._parse_apply_expr()

    def _parse_apply_expr(self):
        # TODO
        return self._parse_get_expr()

    def _parse_get_expr(self):
        # TODO
        return self._parse_term()

    def _parse_term(self):
        t = self.queue[0].type

        if t == TokenType.str:
            return node.String(self._eat().value)

        elif t == TokenType.num:
            value = self._eat().value
            if value.find('.') > 0:
                return node.Number(float(value))
            else:
                return node.Number(int(value))

        elif t == TokenType.id:
            return self._parse_id()

        elif t == TokenType.lbrace:
            return self._parse_record()

        elif t == TokenType.lparen:
            return self._parse_group()

        else:
            raise ValueError(f'Unrecognized token: {self.queue[0]}')

    _parse_simple_expr = _parse_pipe_expr


