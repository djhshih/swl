import os
import unittest as ut

from swl.ir.force import force_file


_ALIGN = '''# @ Align
# in
#   fastq1, fastq2 file
#   ref file
#   ref_fai file
#   outbase str
# out
#   bam file = ${outbase}.bam
echo align
'''

_SORT = '''# @ Sort
# in
#   bam file
#   outbase str
# out
#   bam file = ${outbase}.bam
#   bai file = ${outbase}.bai
echo sort
'''

_PIPE = '''align = import "align.sh"
sort = import "sort.sh"
\\x ->
    a = align x
    sort (x // a)
'''

_PARTIAL = '''align = import "align.sh"
\\x ->
    align { ref: x.ref }
'''

_CHAIN = '''align = import "align.sh"
sort = import "sort.sh"
align | sort
'''


class TestForce(ut.TestCase):
    def _files(self):
        root = os.path.abspath('/virtual-force')
        return {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sort.sh'): _SORT,
            os.path.join(root, 'pipe.swl'): _PIPE,
            os.path.join(root, 'partial.swl'): _PARTIAL,
            os.path.join(root, 'chain.swl'): _CHAIN,
        }, root

    def test_force_saturated_workflow_produces_task_dag(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['name'] for task in data['tasks']], ['align', 'sort'])
        bam = data['tasks'][1]['inputs']['bam']
        self.assertEqual(bam['kind'], 'task_output')
        self.assertEqual(bam['task'], 't1')
        self.assertIn('fastq1', data['inputs'])
        self.assertIn('outbase', data['inputs'])
        self.assertEqual(data['tasks'][0]['task']['doc'], 'Align')
        self.assertIn('echo align', data['tasks'][0]['task']['body'])
        self.assertIn('inputs', data['tasks'][0]['task'])
        self.assertIn('outputs', data['tasks'][0]['task'])
        self.assertEqual(data['tasks'][0]['task']['outputs']['bam']['default']['kind'], 'word')
        self.assertEqual(data['tasks'][0]['deps'], [])
        self.assertEqual(data['tasks'][1]['deps'], ['t1'])

    def test_force_partial_workflow_keeps_function_output(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        data = dag.to_dict()
        self.assertEqual(data['tasks'], [])
        self.assertEqual(data['outputs']['result']['kind'], 'function')

    def test_chain_root_is_instantiated_during_force(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'chain.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['name'] for task in data['tasks']], ['align', 'sort'])
        self.assertIn('fastq1', data['inputs'])
        self.assertEqual(data['tasks'][1]['deps'], ['t1'])

    def test_serialized_dag_is_self_contained(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        text = repr(data)
        self.assertNotIn('Import(', text)
        self.assertNotIn('workflow', text)
        self.assertEqual(data['tasks'][0]['path'], os.path.join(root, 'align.sh'))
        self.assertIn('inputs', data['tasks'][0]['task'])
        self.assertIn('outputs', data['tasks'][0]['task'])


if __name__ == '__main__':
    ut.main()
