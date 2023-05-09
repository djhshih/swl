import ast
import lexer as lex
import queue

class Parser:

    def parse(self, s):
        '''Parse string `s` into an abstract syntax tree.'''
        self.lexer = lex.Lexer(s)
        self.token = next(lexer)
        self.queue = []
        return self.parse_block()

    def eof():
        '''Return whether parser is at end of file.'''
        return self.token.type == lex.TokenType.eof

    def at(self, i = 0):
        '''Peek token at position `i`.'''
        if i == 0:
            # return current token
            return self.token
        elif i < len(self.queue):
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
            self.token = self.queue.pop(0)
        else:
            # otherwise, eat next token
            self.token = next(self.lexer)
        return self.token

    def until(self, token_type):
        '''Return whether current token is of type `token_type`.'''
        return not self.eof() and self.token.type != token_type

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
        while not self.eof():
            exprs.append(self.parse_expr())
        return exprs

    def parse_expr(self):
        t = self.token.type
        if t == lex.TokenType.id:
            return parse_binding()
        else:
            return self.parse_simple_expr()

    # order of precedence for operators
    # from outer-most parser (lowest precedence)
    # to inner-most parser (highest precedence)
    # pipe, update, apply, get

    def parse_pipe_expr():
        # TODO
        return parse_update_expr()

    def parse_update_expr():
        # TODO
        return parse_apply_expr()

    def parse_apply_expr():
        # TODO
        return parse_get_expr()

    def parse_get_expr():
        # TODO
        return ast.ExprN()

    parse_simple_expr = parse_pipe_expr


