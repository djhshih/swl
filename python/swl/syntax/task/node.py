from enum import Enum


class SectionType(Enum):
    IN = 'in'
    OUT = 'out'
    RUN = 'run'


class Task:
    def __init__(self, annotation, body: str):
        self.annotation = annotation
        self.body = body


class Annotation:
    def __init__(self, doc, sections):
        self.doc = doc
        self.sections = sections


class Section:
    def __init__(self, kind: SectionType, params):
        self.kind = kind
        self.params = params


class Param:
    def __init__(self, names, param_type=None, default=None, desc=None):
        self.names = names
        self.type = param_type
        self.default = default
        self.desc = desc
