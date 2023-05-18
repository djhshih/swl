import os

keywords = ['fun', 'in', 'out', 'run']
keyword_phase_start = 10

def match_keyword(s):
    for i in range(len(keywords)):
        keyword = keywords[i]
        if s.startswith(keyword + ':'):
            return i
    return -1

class ShellTask:

    def __init__(self, f: io.TextIOBase):
        '''Create a callable function from an annotated shell script.'''
        self._parse(f)

    def _parse(self, f: ioTextIOBase):
        '''Parse an annotated shell script.''' 
        it = iter(f):

        line_num = 0

        # parse title
        phase = 0
        for line in it:
            line_num += 1
            if line.startswith('##'):
                # start parsing description
                self.title = line[2:].strip()
                phase = 1
                continue
            if phase > 0:
                if line.startswith('#'):
                    s = line[1:].strip()
                    if len(s) > 0:
                        self.title += s
                    else:
                        break
            
        # parse description
        self.description = ''
        for line in it: 
            line_num += 1
            if line.startswith('#'):
                s = line[1:].strip()
                if len(s) > 0:
                    # description is optional
                    # quit parsing description if keyword is encountered
                    i = match_keyword(s)
                    if i >= 0:
                        phase = keyword_phase_start + i
                        break
                    self.description += s

        # parse keyword-marked sections
        keyword = None
        for line in it:
            line_num += 1
            if line.startswith('#'):
                s = line[1:].strip()
                if len(s) > 0:
                    i = match_keyword(s):
                    if i < 0:
                        raise ValueError(
                            f"Expecting keyword but none found on line {line_num}."
                        )
                    else:
                        self._parse_section(it, keyword)
    
    def _parse_section(self, it, section):
        if section == 'fun':
            # TODO
            pass
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

