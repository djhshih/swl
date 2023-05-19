import os

keywords = ['fun', 'in', 'out', 'run']
keyword_phase_start = 10

def match_keyword(s):
    '''Return the matched keyword or None if no match.'''
    for i in range(len(keywords)):
        keyword = keywords[i]
        if s.strip() == keyword + ':'
            return keyword
    return None

class EndOfAnnotation(Exception):
    '''Raised when end of annotation is encountered'''
    pass

class Parameter:
    '''Parameter for a task.'''
    def __init__(self, name, _type = 'file', doc = ''):
        self.name = name
        self.type = _type
        self.doc = doc

class ShellTask:

    def __init__(self, f: io.TextIOBase):
        '''Create a callable function from an annotated shell script.'''
        self.iter = iter(f)
        self.i = 0
        self._parse()

    def _size(self) -> int:
        '''Remaining size of input string'''
        return len(self.lines) - self.i

    def _next_line(self, prefix=None) -> str:
        while True:
            self.i += 1
            line = next(self.iter)
            if prefix:
                if line.startswith(prefix):
                    s = line[1:].strip()
                    if len(s) > 0:
                        return s
                else:
                    raise EndAnnotation
            else:
                if len(line) > 0:
                    return line

    def _parse(self):
        '''Parse an annotated shell script.''' 

        try:
            # parse title
            line = _next_line('##')
            self.title = line[2:].strip()
            
            # parse description,
            # which includes all lines until the next section
            self.description = ''
            k = 0
            while True:
                line = _next_line('#')
                section = match_keyword(line)
                if section:
                    break
                else:
                    self.description += line

            # parse each section
            while True:
                self._parse_section(keywords[k])
                line = _next_line('#')
                section = match_keyword(line)
                if not section:
                    raise ValueError(
                        f"Expecting keyword but none found on line {self.i}."
                    )

        except EndAnnotation:
            pass

        except StopIteration:
            pass

        # unpack the remaining lines
        self.body = [*self.iter]
    
    def _parse_section(self, section):
        line = _next_line('#')

        next_section = match_keyword(line)
        if next_section:
            self._parse_section(next_section)
            return

        if section == 'fun':
            try:
                line = line.strip()
                parts = line.split('->')
                self.inputs  = [Parameter(x.strip()) for x in parts[0].split(',')]
                self.outputs = [Parameter(x.strip()) for x in parts[1].split(',')]
            except:
                raise ValueError(f"fun section is malformed:\n{line}")
        elif section == 'in':
            # TODO
            pass
        elif section == 'out':
            # TODO
            pass
        elif section == 'run':
            # TODO
            pass
        else:
            raise ValueError(f"Unrecognized section: {section}")

