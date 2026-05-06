class Word:
    def __init__(self, parts):
        self.parts = parts

    def __repr__(self):
        return f'Word({self.parts!r})'

    def __eq__(self, other):
        return type(self) is type(other) and self.parts == other.parts


class Literal:
    def __init__(self, text: str):
        self.text = text

    def __repr__(self):
        return f'Literal({self.text!r})'

    def __eq__(self, other):
        return type(self) is type(other) and self.text == other.text


class Var:
    def __init__(self, name: str):
        self.name = name

    def __repr__(self):
        return f'Var({self.name!r})'

    def __eq__(self, other):
        return type(self) is type(other) and self.name == other.name


class Expr:
    def __init__(self, text: str):
        self.text = text

    def __repr__(self):
        return f'Expr({self.text!r})'

    def __eq__(self, other):
        return type(self) is type(other) and self.text == other.text


class Parser:
    def parse_word(self, s: str) -> Word:
        parts = []
        i = 0
        literal = []

        while i < len(s):
            c = s[i]
            if c != '$':
                literal.append(c)
                i += 1
                continue

            if literal:
                parts.append(Literal(''.join(literal)))
                literal = []

            if i + 1 >= len(s):
                literal.append('$')
                i += 1
                continue

            if s[i + 1] == '{':
                end = s.find('}', i + 2)
                if end == -1:
                    raise ValueError('Unterminated interpolation expression')
                text = s[i + 2:end].strip()
                if not text:
                    raise ValueError('Empty interpolation expression')
                if self._is_name(text):
                    parts.append(Var(text))
                else:
                    parts.append(Expr(text))
                i = end + 1
                continue

            j = i + 1
            while j < len(s) and (s[j].isalnum() or s[j] == '_'):
                j += 1
            if j == i + 1:
                literal.append('$')
                i += 1
                continue
            parts.append(Var(s[i + 1:j]))
            i = j

        if literal:
            parts.append(Literal(''.join(literal)))

        return Word(parts)

    def _is_name(self, s: str) -> bool:
        if not s:
            return False
        if not (s[0].isalpha() or s[0] == '_'):
            return False
        for c in s[1:]:
            if not (c.isalnum() or c == '_'):
                return False
        return True


def parse_word(s: str) -> Word:
    return Parser().parse_word(s)
