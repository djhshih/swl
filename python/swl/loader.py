import os


class Loader:
    def __init__(self, files=None):
        self.files = files or {}
        self._parsed_tasks = {}
        self._checked_workflows = {}
        self._loading = []

    def read_file(self, path: str) -> str:
        if path in self.files:
            return self.files[path]
        with open(path, 'r') as f:
            content = f.read()
        self.files[path] = content
        return content

    def get_parsed_task(self, path: str):
        return self._parsed_tasks.get(path)

    def cache_task(self, path: str, task, signature, parsed_body):
        self._parsed_tasks[path] = (task, signature, parsed_body)

    def get_checked_workflow(self, path: str):
        return self._checked_workflows.get(path)

    def cache_workflow(self, path: str, check):
        self._checked_workflows[path] = check
