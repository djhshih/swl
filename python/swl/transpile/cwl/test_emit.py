import os
import unittest

from swl.ir.dag import DAG, Input, Literal, Merge, Record, StepCall
from swl.ir.force import force_file
from swl.transpile.cwl.emit import transpile_dag_dict


class TestCWLTranspile(unittest.TestCase):
    def _files(self):
        root = '/tmp/swl-cwl'
        return {
            os.path.join(root, 'align.sh'): '''#@  Align paired-end sequencing reads
#
# in
#   fastq1, fastq2 file | paired-end reads
#   ref file | reference sequence
#   ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa file
#     | reference bwa index files
#   outbase str | output base name
#
# out
#   bam file = ${outbase}.bam
#     | output alignment
#
# run
#   cpu = 2
#   time = 2-04:30:00
#   image = djhshih/seqkit:0.1

bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam
''',
            os.path.join(root, 'sort.sh'): '''#@  Sort alignment by coordinates and index
# in
#   bam      file                    | input bam
#   outbase  str                     | output base name
# out
#   bam      file  =  ${outbase}.bam | output alignment
#   bai      file  =  ${outbase}.bai | output alignment index
# run
#   cpu    = 2
#   memory = 8G
#   image  = djhshih/seqkit:0.1

samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam
samtools index ${outbase}.bam ${outbase}.bai
''',
            os.path.join(root, 'call.sh'): '''#@  Call mutations on read alignment
# in
#   bam      file                    | input bam
#   ref      file                    | reference sequence
#   ref_fai  file                    | reference index
#   outbase  str                     | output base name
# out
#   bcf      file  =  ${outbase}.bcf | output alignment
# run
#   time   = 30
#   image  = djhshih/seqkit:0.1

bcftools mpileup -Ou -f ${ref} ${bam} | 
	bcftools call -mv -Ob -o ${outbase}.bcf
''',
            os.path.join(root, 'function.swl'): '''align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

\\x ->
    a = align x
    s = sort ( x // a )
    c = call ( x // a // s )
    a // s // c
''',
            os.path.join(root, 'partial.swl'): '''align = import "align.sh"
align_hg38 = align {
    ref: "hg38.fa",
    ref_amb: "hg38.fa.amb",
    ref_ann: "hg38.fa.ann",
    ref_bwt: "hg38.fa.bwt",
    ref_pac: "hg38.fa.pac",
    ref_sa: "hg38.fa.sa",
    cpu: 4
}
align_hg38
''',
            os.path.join(root, 'literal_output.swl'): '1\n',
            os.path.join(root, 'import_partial.swl'): 'partial = import "partial.swl"\npartial\n',
            os.path.join(root, 'bad_expr.sh'): '''#@  Bad expr output
# in
#   outbase str
# out
#   bam file = ${outbase / 2}.bam

echo hi
''',
            os.path.join(root, 'bad_expr.swl'): 'bad = import "bad_expr.sh"\nbad\n',
            os.path.join(root, 'merge.sh'): '''#@  Merge bam
# in
#   bam [file]
#   outbase str
# out
#   bam file = ${outbase}.bam

echo merge
''',
            os.path.join(root, 'batch.swl'): '''align = import "align.sh"
merge = import "merge.sh"

\\xs ->
    calls = map align xs
    merge { bam: calls.bam, outbase: "merged" }
''',
            os.path.join(root, 'mk_align.swl'): '''align = import "align.sh"

\\x ->
    align x
''',
            os.path.join(root, 'batch_workflow.swl'): '''mk = import "mk_align.swl"
merge = import "merge.sh"

\\xs ->
    ys = map mk xs
    merge { bam: ys.bam, outbase: "merged" }
''',
            os.path.join(root, 'batch_lambda.swl'): '''merge = import "merge.sh"

\\xs ->
    f = \\x -> { bam: x.bam }
    ys = map f xs
    merge { bam: ys.bam, outbase: "merged" }
''',
            os.path.join(root, 'sub_lambda.swl'): '''align = import "align.sh"

\\x ->
    a = align x
    { bam2: a.bam }
''',
            os.path.join(root, 'batch_lambda_task.swl'): '''sub = import "sub_lambda.swl"
merge = import "merge.sh"

\\xs ->
    f = \\x -> sub x
    ys = map f xs
    merge { bam: ys.bam2, outbase: "merged" }
''',
        }, root

    def test_transpile_function_workflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        self.assertEqual(cwl['cwlVersion'], 'v1.0')
        self.assertEqual(cwl['$graph'][-1]['id'], '#main')
        tools = [item for item in cwl['$graph'] if item['class'] == 'CommandLineTool']
        self.assertEqual([item['id'] for item in tools], ['#align', '#sort', '#call'])
        workflow = cwl['$graph'][-1]
        outputs = {item['id']: item for item in workflow['outputs']}
        self.assertEqual(outputs['#main/bam']['outputSource'], '#main/sort/bam')
        self.assertEqual(outputs['#main/bai']['outputSource'], '#main/sort/bai')
        self.assertEqual(outputs['#main/bcf']['outputSource'], '#main/call/bcf')
        align = tools[0]
        self.assertEqual(align['baseCommand'], ['bash', 'script.sh'])
        self.assertEqual(align['requirements'][0]['class'], 'InitialWorkDirRequirement')
        self.assertIn('bwa mem', align['requirements'][0]['listing'][0]['entry'])
        self.assertEqual(align['requirements'][1]['coresMin'], 2)
        self.assertEqual(align['requirements'][2]['dockerPull'], 'djhshih/seqkit:0.1')
        self.assertEqual([item['id'] for item in align['inputs']][:2], ['#align/fastq1', '#align/fastq2'])
        self.assertFalse(any('name' in task for task in dag.to_dict()['steps']))

    def test_output_glob_uses_cwl_expression(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        align = next(item for item in cwl['$graph'] if item.get('id') == '#align')
        self.assertEqual(align['outputs'][0]['outputBinding']['glob'], "$(inputs.outbase + '.bam')")

    def test_rejects_non_task_workflow_output(self):
        dag = DAG(inputs={}, steps=[], outputs={'x': Literal(1)})
        with self.assertRaisesRegex(ValueError, 'Unsupported workflow output'):
            transpile_dag_dict(dag.to_dict())

    def test_partial_workflow_transpiles(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        outputs = {item['id']: item for item in workflow['outputs']}
        self.assertEqual(outputs['#main/bam']['outputSource'], '#main/align/bam')

    def test_literal_top_level_workflow_is_rejected_before_transpile(self):
        files, root = self._files()
        with self.assertRaisesRegex(ValueError, 'Workflow must evaluate to a function'):
            force_file(os.path.join(root, 'literal_output.swl'), files)

    def test_imported_workflow_output_transpiles_as_workflow_step(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'import_partial.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        steps = {item['id']: item for item in workflow['steps']}
        self.assertEqual(steps['#main/partial']['run'], '#partial')

    def test_rejects_merged_task_input_binding(self):
        task = StepCall(
            id='align',
            path='/tmp/align.sh',
            bindings={'x': Merge(Input('a'), Input('b'))},
            outputs=['bam'],
            task={
                'body': 'echo hi\n',
                'inputs': {'x': {'type': 'str', 'desc': None}},
                'outputs': {'bam': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'x.bam'}]}, 'desc': None}},
                'run': {},
            },
        )
        bad = {
            'inputs': {'a': {'type': None, 'desc': None}, 'b': {'type': None, 'desc': None}},
            'steps': [
                {
                    'id': 'align',
                    'path': '/tmp/align.sh',
                    'deps': [],
                    'inputs': {'x': {'type': 'str', 'desc': None}},
                    'bindings': {'x': {'source': 'merge', 'left': {'source': 'input', 'name': 'a'}, 'right': {'source': 'input', 'name': 'b'}}},
                    'outputs': {'bam': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'x.bam'}]}, 'desc': None}},
                    'run': {},
                    'script': 'echo hi\n',
                }
            ],
            'outputs': {'bam': {'step': 'align', 'output': 'bam'}},
        }
        with self.assertRaisesRegex(ValueError, 'Unsupported step binding during deserialization: x'):
            transpile_dag_dict(bad)

    def test_batch_mapped_task_emits_scatter_and_array_input_type(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        align = next(step for step in workflow['steps'] if step['id'] == '#main/align')
        self.assertEqual(align['scatterMethod'], 'dotproduct')
        self.assertEqual(sorted(align['scatter']), ['#main/align/fastq1', '#main/align/fastq2', '#main/align/outbase', '#main/align/ref', '#main/align/ref_amb', '#main/align/ref_ann', '#main/align/ref_bwt', '#main/align/ref_pac', '#main/align/ref_sa'])
        merge_tool = next(item for item in cwl['$graph'] if item.get('id') == '#merge')
        bam_input = next(item for item in merge_tool['inputs'] if item['id'] == '#merge/bam')
        self.assertEqual(bam_input['type'], {'type': 'array', 'items': 'File'})

    def test_batch_mapped_workflow_emits_scattered_subworkflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch_workflow.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        mk = next(step for step in workflow['steps'] if step['id'] == '#main/mk')
        self.assertEqual(mk['scatterMethod'], 'dotproduct')
        self.assertEqual(sorted(mk['scatter']), ['#main/mk/fastq1', '#main/mk/fastq2', '#main/mk/outbase', '#main/mk/ref', '#main/mk/ref_amb', '#main/mk/ref_ann', '#main/mk/ref_bwt', '#main/mk/ref_pac', '#main/mk/ref_sa'])
        inputs = {item['id']: item for item in workflow['inputs']}
        self.assertIn('#main/fastq1', inputs)
        self.assertNotIn('#main/xs', inputs)
        subwf = next(item for item in cwl['$graph'] if item.get('id') == '#mk')
        self.assertEqual(subwf['class'], 'Workflow')

    def test_batch_mapped_simple_lambda_emits_generated_scattered_subworkflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch_lambda.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        generated = [item for item in cwl['$graph'] if item.get('class') == 'Workflow' and item.get('id', '').startswith('#map_lambda')]
        self.assertTrue(generated)
        step = next(step for step in workflow['steps'] if step['run'].startswith('#map_lambda'))
        self.assertEqual(step['scatterMethod'], 'dotproduct')

    def test_batch_mapped_lambda_with_inner_task_emits_generated_scattered_subworkflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch_lambda_task.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        generated = next(item for item in cwl['$graph'] if item.get('class') == 'Workflow' and item.get('id', '').startswith('#map_lambda'))
        self.assertTrue(generated['steps'])
        self.assertEqual(generated['steps'][0]['run'], '#sub')
        step = next(step for step in workflow['steps'] if step['run'].startswith('#map_lambda'))
        self.assertEqual(step['scatterMethod'], 'dotproduct')

    def test_rejects_output_expr_interpolation(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'bad_expr.swl'), files)
        with self.assertRaisesRegex(ValueError, 'Unsupported step output path for CWL transpilation: bad.bam: Unsupported interpolation'):
            transpile_dag_dict(dag.to_dict())


if __name__ == '__main__':
    unittest.main()
