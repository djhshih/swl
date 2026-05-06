class Annotation:
    def __init__(self, task_doc, sections):
        self.task_doc = task_doc
        self.sections = sections

class Section:
    def __init__(self, section_type, params):
        self.section_type = section_type  # 'in', 'out', 'run'
        self.params = params

class Parameter:
    def __init__(self, names, param_type=None, default=None, description=None):
        self.names = names  # list of names
        self.type = param_type  # e.g., 'file', 'str', 'int', 'float', 'file?', etc.
        self.default = default  # default value string
        self.description = description  # description string

class AnnotationParser:
    '''Parser for task annotations in bash scripts.'''

    def parse(self, script_content: str) -> Annotation:
        '''Parse annotation comments from a bash script.'''
        lines = script_content.split('\n')
        
        task_doc = None
        sections = []
        current_section = None
        current_params = []
        in_code = False
        
        for line in lines:
            original_line = line
            
            # Check if this is a comment line
            stripped = line.strip()
            if not stripped.startswith('#'):
                # Non-comment line - this is code
                in_code = True
                continue
            
            # Strip leading # and optional space
            if stripped.startswith('#'):
                stripped = stripped[1:]
                if stripped.startswith(' '):
                    stripped = stripped[1:]
            
            # Check for task doc
            if stripped.startswith('@'):
                if current_section is not None and current_params:
                    sections.append(Section(current_section, current_params))
                    current_params = []
                task_doc = stripped[1:].strip()
                continue
            
            # Check for section headers
            if stripped == 'in' or stripped.startswith('in ') or (stripped == '' and current_section != 'in'):
                if stripped.startswith('in ') or stripped == 'in':
                    if current_section is not None and current_params:
                        sections.append(Section(current_section, current_params))
                        current_params = []
                    current_section = 'in'
                    continue
            if stripped == 'out' or stripped.startswith('out '):
                if current_section is not None and current_params:
                    sections.append(Section(current_section, current_params))
                    current_params = []
                current_section = 'out'
                continue
            if stripped == 'run' or stripped.startswith('run '):
                if current_section is not None and current_params:
                    sections.append(Section(current_section, current_params))
                    current_params = []
                current_section = 'run'
                continue
            
            # Parse parameter line (only if we have a current section)
            if current_section and stripped:
                param = self._parse_param(stripped)
                if param:
                    current_params.append(param)
        
        # Add final section
        if current_section is not None and current_params:
            sections.append(Section(current_section, current_params))
        
        return Annotation(task_doc, sections)

    def _parse_param(self, line: str) -> Parameter:
        '''Parse a parameter line: names type? default? desc?'''
        parts = line.split()
        if not parts:
            return None
        
        names = []
        param_type = None
        default = None
        description = None
        
        # Extract description (starts with |)
        desc_idx = -1
        for i, part in enumerate(parts):
            if part.startswith('|'):
                desc_idx = i
                break
        
        if desc_idx >= 0:
            description = ' '.join(parts[desc_idx:])
            parts = parts[:desc_idx]
        
        # Extract default (starts with =)
        default_idx = -1
        for i, part in enumerate(parts):
            if part.startswith('='):
                default_idx = i
                break
        
        if default_idx >= 0:
            default = ''.join(parts[default_idx:])
            parts = parts[:default_idx]
        
        # First part(s) are names (comma-separated)
        # Last part is type (if present)
        if parts:
            # Handle comma-separated names
            name_parts = []
            type_part = None
            
            for part in parts:
                if ',' in part:
                    # Split by comma
                    subparts = part.split(',')
                    for sp in subparts:
                        sp = sp.strip()
                        if sp:
                            name_parts.append(sp)
                    # Check if last subpart has type indicator
                    # This is handled by looking at remaining parts
                else:
                    name_parts.append(part)
            
            # Determine type - it's the last non-name part if present
            # Types: file, str, int, float, [file], etc. with optional ?
            valid_types = {'file', 'str', 'int', 'float', 
                        'file?', 'str?', 'int?', 'float?',
                        '[file]', '[str]', '[int]', '[float]'}
            
            type_candidates = []
            for part in parts:
                for vt in valid_types:
                    if vt in part:
                        type_candidates.append(vt)
            
            if type_candidates:
                # Take the longest match (e.g., [file] over file)
                type_candidates.sort(key=len, reverse=True)
                param_type = type_candidates[0]
                # Remove type from names
                for name in name_parts[:]:
                    if name == param_type or name.endswith(param_type):
                        name_parts.remove(name)
            
            names = name_parts
        
        if not names:
            return None
        
        return Parameter(names, param_type, default, description)


def parse_file(filepath: str) -> Annotation:
    '''Parse annotation from a file.'''
    with open(filepath, 'r') as f:
        content = f.read()
    parser = AnnotationParser()
    return parser.parse(content)