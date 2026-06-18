from dataclasses import dataclass, field

from swl.ir.lower import Lowerer
from swl.semantic.scope import Scope
from swl.dag.evaluator import _assert_canonical, force_value
from swl.dag.finalize import _finalize_dag
from swl.dag.tooldefs import _force_root


@dataclass
class ForceState:
    lowerer: Lowerer
    inputs: dict = field(default_factory=dict)
    steps: list = field(default_factory=list)
    step_cache: dict = field(default_factory=dict)
    tool_defs: dict = field(default_factory=dict)
    call_counter: int = 0
    step_name_counts: dict = field(default_factory=dict)
    variables: dict = field(default_factory=dict)
    forced_variables: dict = field(default_factory=dict)
    apply_cache: dict = field(default_factory=dict)


def make_force_state(files=None):
    return ForceState(lowerer=Lowerer(files=files))


def force(node, files=None, state=None):
    state = state or make_force_state(files=files)
    _assert_canonical(state, node)
    value = force_value(state, node, Scope())
    value = _force_root(state, value)
    return _finalize_dag(state, value)


def force_file(path: str, files=None):
    state = make_force_state(files=files)
    tree = state.lowerer.lower_file(path)
    return force(tree, state=state)
