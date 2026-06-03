from swl.dag.context import ForceEnv
from swl.dag.node import DAG
import sys as _sys


def __getattr__(name):
    if name == 'Forcer':
        from swl.dag.forcer import Forcer as _f
        _sys.modules[__name__].Forcer = _f
        return _f
    if name == 'force_file':
        from swl.dag.forcer import force_file as _f
        _sys.modules[__name__].force_file = _f
        return _f
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
