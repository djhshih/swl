from swl.syntax.task import interpolation


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
