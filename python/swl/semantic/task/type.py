from enum import Enum
from typing import Dict, List, Optional, Set

class TypeKind(Enum):
    FILE = 'file'
    STR = 'str'
    INT = 'int'
    FLOAT = 'float'
    FILE_OPT = 'file?'
    STR_OPT = 'str?'
    INT_OPT = 'int?'
    FLOAT_OPT = 'float?'
    ARRAY_FILE = '[file]'
    ARRAY_STR = '[str]'
    ARRAY_INT = '[int]'
    ARRAY_FLOAT = '[float]'

def parse_type(type_str: str) -> TypeKind:
    '''Parse a type string into a TypeKind.'''
    type_str = type_str.strip()
    
    # Map string to TypeKind
    type_map = {
        'file': TypeKind.FILE,
        'str': TypeKind.STR,
        'int': TypeKind.INT,
        'float': TypeKind.FLOAT,
        'file?': TypeKind.FILE_OPT,
        'str?': TypeKind.STR_OPT,
        'int?': TypeKind.INT_OPT,
        'float?': TypeKind.FLOAT_OPT,
        '[file]': TypeKind.ARRAY_FILE,
        '[str]': TypeKind.ARRAY_STR,
        '[int]': TypeKind.ARRAY_INT,
        '[float]': TypeKind.ARRAY_FLOAT,
    }
    
    if type_str not in type_map:
        raise ValueError(f'Unknown type: {type_str}')
    
    return type_map[type_str]

def is_optional(typ: TypeKind) -> bool:
    '''Check if type is optional.'''
    return typ in (TypeKind.FILE_OPT, TypeKind.STR_OPT, 
                 TypeKind.INT_OPT, TypeKind.FLOAT_OPT)

def is_array(typ: TypeKind) -> bool:
    '''Check if type is an array type.'''
    return typ in (TypeKind.ARRAY_FILE, TypeKind.ARRAY_STR,
                   TypeKind.ARRAY_INT, TypeKind.ARRAY_FLOAT)

def base_type(typ: TypeKind) -> TypeKind:
    '''Get the base type (without optional/array wrapper).'''
    base_map = {
        TypeKind.FILE_OPT: TypeKind.FILE,
        TypeKind.STR_OPT: TypeKind.STR,
        TypeKind.INT_OPT: TypeKind.INT,
        TypeKind.FLOAT_OPT: TypeKind.FLOAT,
        TypeKind.ARRAY_FILE: TypeKind.FILE,
        TypeKind.ARRAY_STR: TypeKind.STR,
        TypeKind.ARRAY_INT: TypeKind.INT,
        TypeKind.ARRAY_FLOAT: TypeKind.FLOAT,
    }
    return base_map.get(typ, typ)

# Type compatibility matrix from spec
COMPATIBLE_PAIRS = {
    # output -> input allowed pairs
    (TypeKind.FILE, TypeKind.FILE): True,
    (TypeKind.FILE, TypeKind.FILE_OPT): True,
    (TypeKind.STR, TypeKind.STR): True,
    (TypeKind.STR, TypeKind.STR_OPT): True,
    (TypeKind.INT, TypeKind.INT): True,
    (TypeKind.INT, TypeKind.INT_OPT): True,
    (TypeKind.FLOAT, TypeKind.FLOAT): True,
    (TypeKind.FLOAT, TypeKind.FLOAT_OPT): True,
}

def types_compatible(output_type: TypeKind, input_type: TypeKind) -> bool:
    '''Check if output type is compatible with input type.'''
    # Direct match
    if (output_type, input_type) in COMPATIBLE_PAIRS:
        return True
    
    # Allow output to satisfy optional input if base types match
    output_base = base_type(output_type)
    input_base = base_type(input_type)
    
    if is_optional(input_type) and output_base == input_base:
        return True
    
    return False


class TaskSignature:
    '''Type signature for a task.'''
    
    def __init__(self, inputs: Dict[str, TypeKind], 
                 outputs: Dict[str, TypeKind],
                 run: Dict[str, any] = None):
        self.inputs = inputs  # name -> TypeKind
        self.outputs = outputs  # name -> TypeKind
        self.run = run or {}  # runtime parameters like cpu, memory, image
    
    def input_names(self) -> Set[str]:
        return set(self.inputs.keys())
    
    def output_names(self) -> Set[str]:
        return set(self.outputs.keys())


class TypeChecker:
    '''Type checker for workflows.'''
    
    def __init__(self):
        self.signatures: Dict[str, TaskSignature] = {}
    
    def add_task(self, name: str, signature: TaskSignature):
        '''Add a task signature.'''
        self.signatures[name] = signature
    
    def check_chain(self, left: str, right: str) -> List[str]:
        '''Check type compatibility when chaining left | right.
        Returns list of error messages (empty if compatible).'''
        errors = []
        
        if left not in self.signatures:
            errors.append(f'Unknown task: {left}')
            return errors
        
        if right not in self.signatures:
            errors.append(f'Unknown task: {right}')
            return errors
        
        left_sig = self.signatures[left]
        right_sig = self.signatures[right]
        
        # Check each output of left against inputs of right
        for out_name, out_type in left_sig.outputs.items():
            if out_name in right_sig.inputs:
                in_type = right_sig.inputs[out_name]
                if not types_compatible(out_type, in_type):
                    errors.append(
                        f'Type mismatch for "{out_name}": '
                        f'{out_type.value} -> {in_type.value}'
                    )
        
        return errors
    
    def infer_workflow_input(self, tasks: List[str]) -> Set[str]:
        '''Infer the required workflow input from a list of tasks.'''
        # Collect all outputs from all tasks
        all_outputs = set()
        for task in tasks:
            if task in self.signatures:
                all_outputs.update(self.signatures[task].output_names())
        
        # Inputs that are not satisfied by any upstream task
        unsatisfied = set()
        for task in tasks:
            if task in self.signatures:
                inputs = self.signatures[task].input_names()
                for inp in inputs:
                    if inp not in all_outputs:
                        unsatisfied.add(inp)
        
        return unsatisfied