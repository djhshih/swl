import json
import os
import tempfile
import unittest

from swl.api import force_workflow, load_workflow, transpile_dag
from swl.dag.node import DAG


class TestAPI(unittest.TestCase):
    def _files(self):
        return {
            '/v/align.sh': (
                '# @ Align\n# in\n#   x file\n# out\n'
                '#   y file = out.txt\necho hi\n'
            ),
            '/v/wf.swl': 'a = import "align.sh"\na\n',
        }

    def test_force_workflow_returns_dag(self):
        dag = force_workflow('/v/wf.swl', files=self._files())
        self.assertIsInstance(dag, DAG)
        self.assertEqual(len(dag.steps), 1)

    def test_load_workflow_returns_check_result(self):
        result = load_workflow('/v/wf.swl', files=self._files())
        self.assertIsNotNone(result.signature)
        self.assertIn('y', result.signature.outputs)

    def test_transpile_dag_cwl(self):
        dag = force_workflow('/v/wf.swl', files=self._files())
        td = tempfile.TemporaryDirectory()
        dag.write(os.path.join(td.name, 'plan.json'))
        result = transpile_dag(os.path.join(td.name, 'plan.json'), 'cwl')
        parsed = json.loads(result)
        self.assertEqual(parsed['cwlVersion'], 'v1.2')

    def test_transpile_dag_nf(self):
        dag = force_workflow('/v/wf.swl', files=self._files())
        td = tempfile.TemporaryDirectory()
        dag.write(os.path.join(td.name, 'plan.json'))
        result = transpile_dag(os.path.join(td.name, 'plan.json'), 'nf')
        self.assertIn('process', result)

    def test_transpile_dag_invalid_target(self):
        with self.assertRaises(ValueError):
            transpile_dag('/fake.json', 'invalid')


if __name__ == '__main__':
    unittest.main()
