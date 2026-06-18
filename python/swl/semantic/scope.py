class DuplicateBindingError(Exception):
    pass


class Binding:
    def __init__(self, name, node=None, value=None, **kwargs):
        self.name = name
        self.node = node
        self.value = value
        for key, value in kwargs.items():
            setattr(self, key, value)


class Scope:
    def __init__(self, parent=None):
        self.parent = parent
        self.locals = {}

    def declare(self, name, node=None, value=None, **kwargs):
        if name in self.locals:
            raise DuplicateBindingError(name)
        binding = Binding(name, node=node, value=value, **kwargs)
        self.locals[name] = binding
        return binding

    def set_local(self, name, value=None, **kwargs):
        if name in self.locals:
            binding = self.locals[name]
            if value is not None:
                binding.value = value
            for key, val in kwargs.items():
                setattr(binding, key, val)
            return binding
        return self.declare(name, value=value, **kwargs)

    def resolve(self, name):
        if name in self.locals:
            return self.locals[name]
        if self.parent is not None:
            return self.parent.resolve(name)
        return None

    def __contains__(self, name):
        return self.resolve(name) is not None

    def __or__(self, other):
        if isinstance(other, set):
            result = Scope(parent=self.parent)
            result.locals.update(self.locals)
            for item in other:
                if item not in result.locals:
                    result.declare(item)
            return result
        if isinstance(other, Scope):
            result = Scope(parent=self.parent)
            result.locals.update(self.locals)
            result.locals.update(other.locals)
            return result
        return NotImplemented
