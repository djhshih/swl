_SENTINEL = object()


class ForceEnv:
    def __init__(self, parent=None, values=None):
        self.parent = parent
        self.values = values or {}

    def bind(self, name, value):
        self.values[name] = value

    def lookup(self, name):
        if name in self.values:
            return self.values[name]
        if self.parent is not None:
            return self.parent.lookup(name)
        return None
