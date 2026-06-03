import os

from swl.semantic.task.type import TypeChecker
from swl.semantic.wf.imports import load_import, read_file, load_imports
from swl.semantic.wf.infer import infer_inputs
from swl.semantic.wf.scope import check_scope, check_chains
from swl.semantic.wf.signature import build_workflow_signature
from swl.syntax.wf.parser import Parser as WfParser


class WorkflowCheck:
    def __init__(self, tree, imports, errors, inferred_inputs, signature=None, is_batch=False, workflow_type=None, root_input_type=None, root_output_type=None):
        self.tree = tree
        self.imports = imports
        self.errors = list(errors)
        self.inferred_inputs = inferred_inputs
        self.signature = signature
        self.is_batch = is_batch
        self.workflow_type = workflow_type
        self.root_input_type = root_input_type
        self.root_output_type = root_output_type


class Checker:
    def __init__(self, files=None):
        self._loading = []
        self.files = files or {}
        self._apply_satisfied = {}

    def _load_import(self, name: str, path: str):
        return load_import(self, name, path)

    def load(self, path: str) -> WorkflowCheck:
        full_path = os.path.abspath(path)
        if full_path in self._loading:
            raise ValueError(f'Circular workflow import: {full_path}')
        self._loading.append(full_path)
        try:
            src = read_file(self, full_path)
            tree = WfParser().parse(src)
            imports = load_imports(self, tree, os.path.dirname(full_path))
            errors = check_scope(self, tree)
            checker = TypeChecker()
            for imported in imports.values():
                checker.add_task(imported.name, imported.signature)
            self._type_checker = checker
            errors.extend(check_chains(self, tree, checker))
            inferred_inputs, infer_errors = infer_inputs(self, tree, imports)
            signature, is_batch, workflow_type, root_input_type, root_output_type = build_workflow_signature(self, tree, imports, inferred_inputs, errors)
            errors.extend(infer_errors)
            if signature is None:
                errors.append('Workflow must evaluate to a function')
            return WorkflowCheck(tree, imports, errors, inferred_inputs, signature, is_batch, workflow_type, root_input_type, root_output_type)
        finally:
            self._loading.pop()

    def load_content(self, content: str, path: str = '<memory>.swl') -> WorkflowCheck:
        full_path = os.path.abspath(path)
        old = self.files.get(full_path)
        self.files[full_path] = content
        try:
            return self.load(full_path)
        finally:
            if old is None:
                del self.files[full_path]
            else:
                self.files[full_path] = old
