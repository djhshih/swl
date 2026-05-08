import os
import unittest as ut

from swl.ir.force import force_file


_ALIGN = '''# @ Align
# in
#   fastq1, fastq2 file
#   | fastq reads
#   ref file
#   | reference fasta
#   ref_fai file
#   | reference index
#   outbase str
#   | output base
# out
#   bam file = ${outbase}.bam
#   | aligned bam
echo align
'''

_SORT = '''# @ Sort
# in
#   bam file
#   | input bam
#   outbase str
#   | output base
# out
#   bam file = ${outbase}.bam
#   | sorted bam
#   bai file = ${outbase}.bai
#   | sorted index
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

_CALL = '''# @ Call
# in
#   bam file
#   | input bam
#   ref file
#   | reference sequence
#   ref_fai file
#   | reference index
#   outbase str
#   | output base
# out
#   bcf file = ${outbase}.bcf
#   | called variants
echo call
'''

_CHAIN = '''align = import "align.sh"
sort = import "sort.sh"
call = import "call.sh"
align | sort | call
'''


class TestForce(ut.TestCase):
    def _files(self):
        root = os.path.abspath('/virtual-force')
        return {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sort.sh'): _SORT,
            os.path.join(root, 'call.sh'): _CALL,
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
        self.assertEqual(bam['source'], 'task')
        self.assertEqual(bam['task'], 't1')
        self.assertIn('fastq1', data['inputs'])
        self.assertIn('outbase', data['inputs'])
        self.assertEqual(data['inputs']['fastq1']['type'], 'file')
        self.assertEqual(data['inputs']['fastq1']['desc'], 'fastq reads')
        self.assertEqual(data['inputs']['outbase']['type'], 'str')
        self.assertEqual(data['inputs']['outbase']['desc'], 'output base')
        self.assertIn('echo align', data['tasks'][0]['script'])
        self.assertIn('outputs', data['tasks'][0])
        self.assertEqual(data['tasks'][0]['outputs']['bam']['default']['kind'], 'word')
        self.assertEqual(data['tasks'][0]['deps'], [])
        self.assertEqual(data['tasks'][1]['deps'], ['t1'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bai', 'bam'])

    def test_force_partial_workflow_keeps_function_output(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        data = dag.to_dict()
        self.assertEqual(data['tasks'], [])
        self.assertEqual(data['outputs']['result']['source'], 'function')

    def test_chain_root_is_instantiated_during_force(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'chain.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['name'] for task in data['tasks']], ['align', 'sort', 'call'])
        self.assertIn('fastq1', data['inputs'])
        self.assertIn('ref_fai', data['inputs'])
        self.assertEqual(data['inputs']['ref_fai']['type'], 'file')
        self.assertEqual(data['tasks'][2]['deps'], ['t2'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bai', 'bam', 'bcf'])
        self.assertEqual(data['outputs']['bam']['source'], 'task')
        self.assertEqual(data['outputs']['bai']['source'], 'task')
        self.assertEqual(data['outputs']['bcf']['source'], 'task')

    def test_serialized_dag_is_self_contained(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        text = repr(data)
        self.assertNotIn('Import(', text)
        self.assertNotIn('workflow', text)
        self.assertEqual(data['tasks'][0]['path'], os.path.join(root, 'align.sh'))
        self.assertIn('script', data['tasks'][0])
        self.assertIn('outputs', data['tasks'][0])


if __name__ == '__main__':
    ut.main()
