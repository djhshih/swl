from swl.parse import node
from swl.parse import lexer as lex
from swl.parse.lexer import TokenType

class Parser:

    def parse(self, s):
        '''Parse string `s` into an abstract syntax tree.'''
        self.lexer = lex.Lexer(s)
        self.queue = [next(self.lexer)]
        return self.parse_block()

    def eof(self):
        '''Return whether parser is at end of file.'''
        if not self.queue: return True
        return self.queue[0].type == TokenType.eof

    def at(self, i = 0):
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
    
    def eat(self):
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

    def expect(self, token_type):
        '''Eat a token and check that it has the expected type'''
        token = self.eat()
        if not token:
            raise ValueError(f'Parsing failed due to undefined token')
        if token.type != token_type:
            raise ValueError(
                f'Parsing failed due to unexpected token: expected token type {token_type.name}, but got token {token}'
            )
        return token

    def until(self, token_type):
        '''Return whether current token is of type `token_type`.'''
        return not self.eof() and self.qu.type != token_type

    def find(self, token_type):
        '''Find and return index of token with type `token_type`,
           but do not advance.'''
        i = 0
        while not self.eof():
            token = self.at(i)
            if token.type == token_type:
                return i
            i += 1
        return -1

    def parse_block(self):
        exprs = []
        # TODO or check for dedent
        while True:
            exprs.append(self.parse_expr())
            if self.eof(): break
            self.expect(TokenType.eol)
        # TODO eat the dedent
        return node.Block(exprs)

    def parse_expr(self):
        t = self.queue[0].type

        if t == TokenType.bslash:
            return self.parse_function()

        elif t == TokenType.id and self.at(1).type == TokenType.equal:
            return self.parse_binding()

        else:
            return self.parse_simple_expr()

    # order of precedence for operators
    # from outer-most parser (lowest precedence)
    # to inner-most parser (highest precedence)
    # pipe, update, apply, get

    def parse_id(self):
        iden = self.expect(TokenType.id)
        return node.Identifier(iden.value)

    def parse_record(self):
        # TODO
        return node.Record({})

    def parse_group(self):
        self.eat()  # eat open parenthesis
        value = self.parse_expr()
        self.expect(TokenType.rparen)
        return value

    def parse_function(self):
        self.eat()  # eat backslash
        param = self.parse_id()
        self.expect(TokenType.arrow)
        block = self.parse_block()
        return node.Function(iden, block)

    def parse_binding(self):
        print('queue:', self.queue)
        iden = self.parse_id()
        self.expect(TokenType.equal)
        expr = self.parse_expr()
        return node.Binding(iden, expr)

    def parse_pipe_expr(self):
        # TODO
        return self.parse_update_expr()

    def parse_update_expr(self):
        # TODO
        return self.parse_apply_expr()

    def parse_apply_expr(self):
        # TODO
        return self.parse_get_expr()

    def parse_get_expr(self):
        # TODO
        return self.parse_term()

    def parse_term(self):
        t = self.queue[0].type

        if t == TokenType.str:
            return node.String(self.eat().value)

        elif t == TokenType.num:
            value = self.eat().value
            if value.find('.') > 0:
                return node.Number(float(value))
            else:
                return node.Number(int(value))

        elif t == TokenType.id:
            return self.parse_id()

        elif t == TokenType.lbrace:
            return self.parse_record()

        elif t == TokenType.lparen:
            return self.parse_group()

        else:
            raise ValueError(f'Unexpected token: {self.queue[0]}')

    parse_simple_expr = parse_pipe_expr


