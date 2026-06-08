from swl.syntax.task import node
from swl.syntax.task import interpolation


_VALID_TYPES = {
    'file', 'str', 'int', 'float',
    'file?', 'str?', 'int?', 'float?',
    '[file]', '[str]', '[int]', '[float]',
}

_SECTION_TYPES = {
    'in': node.SectionType.IN,
    'out': node.SectionType.OUT,
    'run': node.SectionType.RUN,
}


class Parser:
    def parse(self, script: str) -> node.Task:
        annotation_lines, body = self._split_script(script)
        annotation = self._parse_annotation(annotation_lines)
        return node.Task(annotation, body)

    def _split_script(self, script: str):
        lines = script.splitlines()
        annotation = []
        body_start = len(lines)

        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith('#'):
                annotation.append(self._strip_comment(line))
                continue
            if stripped == '':
                if annotation:
                    annotation.append('')
                continue
            body_start = i
            break

        body = '\n'.join(lines[body_start:])
        return annotation, body

    def _strip_comment(self, line: str) -> str:
        stripped = line.lstrip()
        stripped = stripped[1:]
        if stripped.startswith(' '):
            stripped = stripped[1:]
        return stripped.rstrip()

    def _parse_annotation(self, lines):
        self.lines = lines
        self.i = 0

        doc = self._parse_doc()
        sections = []
        while self._skip_blank_lines():
            sections.append(self._parse_section())

        if not sections:
            raise ValueError('Task annotation must contain at least one section')

        return node.Annotation(doc, sections)

    def _parse_doc(self):
        self._skip_leading_blanks()
        if self._eof() or not self._at().startswith('@'):
            raise ValueError('Task annotation must start with a doc line')
        doc = self._eat()[1:].strip()
        return doc

    def _parse_section(self):
        header = self._eat().strip()
        if header not in _SECTION_TYPES:
            raise ValueError(f'Unrecognized section header: {header}')

        params = []
        section_kind = _SECTION_TYPES[header]
        while not self._eof():
            line = self._at().strip()
            if not line:
                self._eat()
                continue
            if line in _SECTION_TYPES:
                break
            if line.startswith('|'):
                if not params:
                    raise ValueError('Description continuation without parameter')
                extra = line[1:].strip()
                if params[-1].desc:
                    params[-1].desc += '\n' + extra
                else:
                    params[-1].desc = extra
                self._eat()
                continue
            params.append(self._parse_param(self._eat(), section_kind))

        return node.Section(section_kind, params)

    def _parse_param(self, line: str, section_kind=None) -> node.Param:
        desc = None
        if '|' in line:
            line, desc = line.split('|', 1)
            desc = desc.strip()

        default = None
        if '=' in line:
            line, default_text = line.split('=', 1)
            default_text = default_text.strip()
            if default_text:
                default = interpolation.parse_word(default_text)

        parts = line.split()
        if not parts:
            raise ValueError('Parameter line is empty')

        param_type = None
        if parts[-1] in _VALID_TYPES:
            param_type = parts.pop()

        names_text = ' '.join(parts).strip()
        names = [x.strip() for x in names_text.split(',') if x.strip()]
        if not names:
            raise ValueError('Parameter line must contain at least one name')

        if param_type is None and section_kind in (node.SectionType.IN, node.SectionType.OUT):
            raise ValueError(
                f'in/out parameter must have a type annotation; '
                f'got {names_text!r} (expected e.g. "str {names_text}" or "file {names_text}")'
            )

        return node.Param(names, param_type, default, desc)

    def _skip_leading_blanks(self):
        while not self._eof() and not self._at().strip():
            self._eat()

    def _skip_blank_lines(self):
        self._skip_leading_blanks()
        return not self._eof()

    def _eof(self):
        return self.i >= len(self.lines)

    def _at(self):
        return self.lines[self.i]

    def _eat(self):
        line = self.lines[self.i]
        self.i += 1
        return line


def parse_file(path: str) -> node.Task:
    with open(path, 'r') as f:
        return Parser().parse(f.read())
