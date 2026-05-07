import os
import unittest as ut

from swl.semantic.wf.check import Checker


_TESTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../../../tests')
)


class TestWorkflowCheck(ut.TestCase):

    def test_load_pipe_workflow(self):
        result = Checker().load(os.path.join(_TESTS_DIR, 'pipe.swl'))
        self.assertEqual(sorted(result.imports.keys()), ['align', 'call', 'sort'])
        self.assertEqual(result.chain_errors, [])
        self.assertEqual(result.inferred_inputs, set())

    def test_load_function_workflow_infers_inputs(self):
        result = Checker().load(os.path.join(_TESTS_DIR, 'function.swl'))
        self.assertEqual(result.chain_errors, [])
        self.assertIn('fastq1', result.inferred_inputs)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)

    def test_load_explicit_workflow_infers_inputs(self):
        result = Checker().load(os.path.join(_TESTS_DIR, 'explicit.swl'))
        self.assertEqual(result.chain_errors, [])
        self.assertIn('fastq1', result.inferred_inputs)
        self.assertIn('fastq2', result.inferred_inputs)
        self.assertIn('ref', result.inferred_inputs)
        self.assertIn('ref_fai', result.inferred_inputs)
        self.assertIn('outbase', result.inferred_inputs)

    def test_bad_pipe_reports_chain_error(self):
        result = Checker().load(os.path.join(_TESTS_DIR, 'bad_pipe.swl'))
        self.assertEqual(len(result.chain_errors), 1)
        self.assertIn('outbase', result.chain_errors[0])


if __name__ == '__main__':
    ut.main()
