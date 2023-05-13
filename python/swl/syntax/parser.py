from swl.syntax import node
from swl.syntax.node import NodeType
from swl.syntax import lexer as lex
from swl.syntax.lexer import Token, TokenType

from typing import List

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
                f'Parsing failed due to unexpected token: \
                  expected token type {token_type.name}, but got token {token}'
            )
        return token

    def _until(self, token_type: TokenType) -> bool:
        '''Return True until eof or current token is of type `token_type`.'''
        return not self.eof() and self.queue[0].type != token_type

    def _until_any(self, token_types: List[TokenType]) -> bool:
        '''Return True until eof or current token type matches any type
           in `token_types`.'''
        if self.eof(): return False
        t = self.queue[0].type
        for token_type in token_types:
            if t == token_type:
                return False
        return True

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

    def _parse_block(self, inner=False) -> node.Expr:
        '''Parse a block of code, which can be outer-level or inner-level.'''
        exprs = []
        if inner:
            # check for block start token
            self._expect(TokenType.bstart)
        while not self.eof():
            exprs.append(self._parse_expr())
            # check for early break because 
            # below tokens are optional right before eof
            if self.eof(): break
            self._expect(TokenType.eol)
            # check for block end token
            if inner and self._at().type == TokenType.bend:
                self._eat()
                break
        return node.Block(exprs)

    def _parse_expr(self) -> node.Expr:
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

    def _parse_id(self) -> node.Expr:
        iden = self._expect(TokenType.id)
        return node.Identifier(iden.value)

    def _parse_record(self) -> node.Expr:
        self._eat()  # eat open brace
        self._optional(TokenType.bstart)
        d = {}
        while self._until_any([TokenType.rbrace, TokenType.bend, TokenType.eol]):
            iden = self._parse_id()
            self._expect(TokenType.colon)
            value = self._parse_expr()
            # comma is required after each key-value pair
            # unless we just parsed the final key-value pair
            if self._at().type != TokenType.rbrace and \
                self._at().type != TokenType.eol:
                self._expect(TokenType.comma)
            d[iden.name] = value
        self._optional(TokenType.eol)
        self._optional(TokenType.bend)
        self._expect(TokenType.rbrace)
        return node.Record(d)

    def _parse_group(self) -> node.Expr:
        self._eat()  # eat open parenthesis
        value = self._parse_expr()
        self._expect(TokenType.rparen)
        return value

    def _parse_function(self) -> node.Expr:
        self._eat()  # eat backslash
        param = self._parse_id()
        self._expect(TokenType.arrow)
        if self._at().type == TokenType.bstart:
            body = self._parse_block(True)
        else:
            body = self._parse_expr()
        return node.Function(param, body)

    def _parse_binding(self) -> node.Expr:
        iden = self._parse_id()
        self._expect(TokenType.equal)
        expr = self._parse_expr()
        return node.Binding(iden, expr)

    def _parse_pipe_expr(self) -> node.Expr:
        left = self._parse_update_expr()

        # only the below node types can potentially have a function
        if \
            left.type != NodeType.fun and \
            left.type != NodeType.id and \
            left.type != NodeType.get and \
            left.type != NodeType.apply and \
            left.type != NodeType.pipe:
            return left

        # only allow pipe if next token is a permissible operand
        # (anything that can evaluate to a function)
        # (a record may have a member function)
        while self._at().type == TokenType.pipe:
            self._eat()  # eat pipe token
            right = self._parse_update_expr()
            left = node.Pipe(left, right)

        return left

    def _parse_update_expr(self) -> node.Expr:
        left = self._parse_apply_expr()

        # only record or identifier can use update notation
        if left.type != NodeType.rec and left.type != NodeType.id:
            return left

        while self._at().type == TokenType.update:
            self._eat()  # eat update token
            right = self._parse_apply_expr()
            left = node.Update(left, right)

        return left

    def _parse_apply_expr(self) -> node.Expr:
        left = self._parse_get_expr()

        # only the below node types can potentially have a callable function
        if \
            left.type != NodeType.fun and \
            left.type != NodeType.id and \
            left.type != NodeType.get and \
            left.type != NodeType.apply and \
            left.type != NodeType.pipe:
            return left

        # only apply function if next token is a permissible operand
        # (anything that can evaluate to a record)
        while \
            self._at().type == TokenType.id or \
            self._at().type == TokenType.lbrace or \
            self._at().type == TokenType.str or \
            self._at().type == TokenType.lparen:
            right = self._parse_get_expr()
            left = node.Apply(left, right)

        return left

    def _parse_get_expr(self) -> node.Expr:
        obj = self._parse_term()

        # only record or identifier can use dot notation
        if obj.type != NodeType.rec and obj.type != NodeType.id:
            return obj

        while self._at().type == TokenType.dot:
            self._eat()  # eat dot token
            member = self._parse_id()
            obj = node.Get(obj, member)

        return obj

    def _parse_term(self) -> node.Expr:
        t = self.queue[0].type

        if t == TokenType.str:
            return node.String(self._eat().value)

        elif t == TokenType.num:
            return node.Number(self._eat().value)

        elif t == TokenType.id:
            return self._parse_id()

        elif t == TokenType.lbrace:
            return self._parse_record()

        elif t == TokenType.lparen:
            return self._parse_group()

        else:
            raise ValueError(f'Unrecognized token: {self.queue[0]}')

    _parse_simple_expr = _parse_pipe_expr


