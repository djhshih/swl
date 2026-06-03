from swl.ir.dag import DAG, Field, ForcedFunction, Input, Literal, Merge, Record, StepCall
from swl.ir.lower import Lowerer, lower_file, parse_and_lower
import sys as _sys


def __getattr__(name):
    if name == 'Forcer':
        from swl.ir.force import Forcer as _f
        _sys.modules[__name__].Forcer = _f
        return _f
    if name == 'force_file':
        from swl.ir.force import force_file as _f
        _sys.modules[__name__].force_file = _f
        return _f
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
