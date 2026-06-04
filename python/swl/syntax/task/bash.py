import re

from swl.syntax.task import interpolation

_BUILTIN_VARS = frozenset({
    'HOME', 'PATH', 'USER', 'SHELL', 'PWD', 'RANDOM',
    'LINENO', 'SECONDS', 'UID', 'HOSTNAME', 'OSTYPE',
    'BASH', 'BASH_VERSION', 'BASH_ENV', 'IFS', 'PS1',
    'PS2', 'PS3', 'PS4', 'OLDPWD', 'SHLVL', 'TERM',
    'LANG', 'LC_ALL', 'LC_CTYPE', 'DISPLAY', 'TMPDIR',
    'EDITOR', 'VISUAL', 'PAGER',
})


def _extract_expr_vars(text: str) -> list[str]:
    names = []
    for m in re.finditer(r'[A-Za-z_]\w*', text):
        names.append(m.group(0))
    return names


def iter_var_refs(script):
    for stmt in script.statements:
        if isinstance(stmt, Command):
            for word in stmt.words:
                for part in _word_parts(word):
                    yield part
        elif isinstance(stmt, Assignment):
            for part in _word_parts(stmt.value):
                yield part


def _word_parts(word):
    parts = word.parts if isinstance(word, interpolation.Word) else [word]
    for part in parts:
        if isinstance(part, interpolation.Var):
            yield (part.name, False)
        elif isinstance(part, interpolation.Expr):
            for name in _extract_expr_vars(part.text):
                yield (name, True)


class Assignment:
    def __init__(self, name: str, value):
        self.name = name
        self.value = value

    def __repr__(self):
        return f'Assignment({self.name!r}, {self.value!r})'

    def __eq__(self, other):
        return type(self) is type(other) and \
            self.name == other.name and self.value == other.value


class Command:
    def __init__(self, text: str, words):
        self.text = text
        self.words = words

    def __repr__(self):
        return f'Command({self.text!r}, {self.words!r})'

    def __eq__(self, other):
        return type(self) is type(other) and \
            self.text == other.text and self.words == other.words


class Script:
    def __init__(self, statements):
        self.statements = statements


class Parser:
    def parse(self, body: str) -> Script:
        statements = []
        lines = body.splitlines()
        while lines and not lines[0].strip():
            lines.pop(0)
        while lines and not lines[-1].strip():
            lines.pop()
        for raw_line in lines:
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue
            assignment = self._parse_assignment(line)
            if assignment is not None:
                statements.append(assignment)
            else:
                statements.append(Command(line, self._parse_words(line)))
        return Script(statements)

    def _parse_assignment(self, line: str):
        if line.startswith('export '):
            line = line[7:].strip()
        parts = line.split(None, 1)
        head = parts[0]
        if '=' not in head:
            return None
        name, value = head.split('=', 1)
        if not self._is_name(name):
            return None
        return Assignment(name, interpolation.parse_word(value))

    def _parse_words(self, line: str):
        words = []
        for part in self._split_shell_words(line):
            if '$' in part:
                words.append(interpolation.parse_word(part))
        return words

    def _split_shell_words(self, line: str):
        words = []
        current = []
        brace_depth = 0

        i = 0
        while i < len(line):
            c = line[i]
            if c.isspace() and brace_depth == 0:
                if current:
                    words.append(''.join(current))
                    current = []
                i += 1
                continue

            if c == '$' and i + 1 < len(line) and line[i + 1] == '{':
                brace_depth += 1
                current.append(c)
                i += 1
                current.append(line[i])
                i += 1
                continue

            if c == '}' and brace_depth > 0:
                brace_depth -= 1

            current.append(c)
            i += 1

        if current:
            words.append(''.join(current))

        return words

    def _is_name(self, s: str) -> bool:
        if not s:
            return False
        if not (s[0].isalpha() or s[0] == '_'):
            return False
        for c in s[1:]:
            if not (c.isalnum() or c == '_'):
                return False
        return True


def parse(body: str) -> Script:
    return Parser().parse(body)
