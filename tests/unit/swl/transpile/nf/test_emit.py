import os
import unittest

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.dag.forcer import force_file
from swl.transpile.nf.emit import transpile_dag_dict


class TestNFTranspile(unittest.TestCase):
    def _files(self):
        root = '/tmp/swl-nf'
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
#   bam      file  =  ${outbase}.sorted.bam | output alignment
#   bai      file  =  ${outbase}.bai | output alignment index
# run
#   cpu    = 2
#   memory = 8G
#   image  = djhshih/seqkit:0.1

samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.sorted.bam
samtools index ${outbase}.sorted.bam ${outbase}.bai
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
            os.path.join(root, 'pipe.swl'): '''align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

align | sort | call
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
            os.path.join(root, 'map_root.swl'): '''call_variant = import "pipe.swl"

map call_variant
''',
        }, root

    def _assert_nf_contains(self, nf, *patterns):
        for pattern in patterns:
            self.assertIn(pattern, nf, f'Expected NF to contain: {pattern}')

    def test_transpile_function_workflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf,
            'process ALIGN {',
            'process SORT {',
            'process CALL {',
            'workflow {',
        )
        self._assert_nf_contains(nf, 'path fastq1')
        self._assert_nf_contains(nf, 'path fastq2')
        self._assert_nf_contains(nf, 'path ref')
        self._assert_nf_contains(nf, 'val outbase')
        self._assert_nf_contains(nf, 'script:')
        self._assert_nf_contains(nf, 'bwa mem -t ${task.cpus}')
        self._assert_nf_contains(nf, '${task.memory / task.cpus}')
        self._assert_nf_contains(nf, 'output:')
        self._assert_nf_contains(nf, 'emit: bam')
        self._assert_nf_contains(nf, 'cpus 2')

    def test_output_interpolation_uses_nf_templates(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self.assertIn('${outbase}.bam', nf)

    def test_rejects_old_style_non_outputspec_workflow_output(self):
        bad = {
            'inputs': {},
            'steps': [],
            'outputs': {'x': {'source': 'literal', 'value': 1}},
        }
        with self.assertRaisesRegex(ValueError, 'has no explicit type'):
            transpile_dag_dict(bad)

    def test_partial_workflow_transpiles(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf, 'workflow {')
        self._assert_nf_contains(nf, 'Channel.fromPath')
        self._assert_nf_contains(nf, 'Channel.value')

    def test_literal_top_level_workflow_is_rejected_before_transpile(self):
        files, root = self._files()
        with self.assertRaisesRegex(ValueError, 'Workflow must evaluate to a function'):
            force_file(os.path.join(root, 'literal_output.swl'), files)

    def test_imported_workflow_inlines_process(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'import_partial.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf, 'ALIGN(')
        self.assertNotIn('MAIN_PARTIAL(', nf)

    def test_rejects_merged_task_input_binding(self):
        bad = {
            'inputs': {'a': {'type': None, 'desc': None}, 'b': {'type': None, 'desc': None}},
            'steps': [{
                'id': 'align',
                'type': 'task',
                'path': '/tmp/align.sh',
                'deps': [],
                'inputs': {'x': {'type': 'str', 'desc': None}},
                'bindings': {'x': {'source': 'merge', 'left': {'source': 'input', 'name': 'a'}, 'right': {'source': 'input', 'name': 'b'}}},
                'outputs': {'bam': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'x.bam'}]}, 'desc': None}},
                'run': {},
                'script': 'echo hi\n',
            }],
            'outputs': {'bam': {'type': 'file', 'desc': None, 'value': {'step': 'align', 'output': 'bam'}}},
        }
        with self.assertRaisesRegex(ValueError, 'Merge'):
            transpile_dag_dict(bad)

    def test_batch_mapped_task_emits_join_and_tuple(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf,
            '.join(',
            'ALIGN(',
            'emit: bam',
        )
        self.assertIn('tuple', nf)
        self.assertIn('path(fastq1)', nf)

    def test_batch_mapped_workflow_inlines_pipeline(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch_workflow.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf, 'ALIGN(fastq1_ch, fastq2_ch,')
        self._assert_nf_contains(nf, 'MERGE(ALIGN.out.bam.toList()')

    def test_batch_mapped_simple_lambda_emits_no_join(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch_lambda.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf, 'MERGE(')
        self.assertNotIn('.join(', nf)

    def test_map_by_transpile_emits_groupTuple(self):
        bad = {
            'inputs': {'sample': {'type': '[str]', 'desc': None}},
            'steps': [{
                'id': 'grouped',
                'type': 'workflow',
                'path': '/tmp/grouped.swl',
                'map': {'source': {'source': 'table', 'name': 'table', 'columns': {'sample': {'source': 'input', 'name': 'sample'}}}, 'group_by': 'sample'},
                'input_schema': {'sample': 'str'},
                'output_schema': {'sample': 'str'},
                'deps': [],
                'inputs': {'sample': {'type': 'str', 'desc': None}},
                'bindings': {},
                'outputs': {'sample': {'type': 'str'}},
                'run': {},
                'script': '',
                'definition': {'class': 'Workflow', 'dag': {'inputs': {'sample': {'type': 'str', 'desc': None}}, 'steps': [], 'outputs': {'sample': {'type': 'str', 'desc': None, 'value': {'source': 'input', 'name': 'sample'}}}}, 'inputs': {'sample': {'type': 'str', 'desc': None}}, 'outputs': {'sample': {'type': 'str'}}, 'body': '', 'run': {}},
            }],
            'outputs': {'sample': {'type': '[str]', 'desc': None, 'value': {'step': 'grouped', 'output': 'sample'}}},
        }
        nf = transpile_dag_dict(bad)
        self._assert_nf_contains(nf, 'groupTuple')

    def test_root_partial_map_inlines_pipeline(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'map_root.swl'), files)
        nf = transpile_dag_dict(dag.to_dict())
        self.assertIn('ALIGN', nf)
        self.assertIn('SORT', nf)
        self.assertIn('CALL', nf)
        self.assertIn('emit: bam', nf)
        self.assertNotIn('.join(', nf)

    def test_process_name_sanitization(self):
        from swl.transpile.nf.emit import _process_name
        self.assertEqual(_process_name('align'), 'ALIGN')
        self.assertEqual(_process_name('Align-BAM'), 'ALIGN_BAM')
        self.assertEqual(_process_name('123align'), '_123ALIGN')
        self.assertEqual(_process_name('_align'), 'ALIGN')

    def test_input_channel_creation(self):
        from swl.transpile.nf.emit import _input_channel
        spec_file = type('Spec', (), {'type': 'file'})()
        spec_str = type('Spec', (), {'type': 'str'})()
        spec_arr = type('Spec', (), {'type': '[file]'})()
        self.assertIn('Channel.fromPath', _input_channel('x', spec_file))
        self.assertIn('Channel.value', _input_channel('x', spec_str))
        self.assertIn('toList', _input_channel('x', spec_arr))

    def test_input_qualifier(self):
        from swl.types import to_nf_qualifier
        self.assertEqual(to_nf_qualifier('file'), ('path', None))
        self.assertEqual(to_nf_qualifier('str'), ('val', 'string'))
        self.assertEqual(to_nf_qualifier('int'), ('val', 'integer'))
        self.assertEqual(to_nf_qualifier('float'), ('val', 'float'))
        self.assertEqual(to_nf_qualifier('[file]'), ('path', None))

    def test_empty_inputs_outputs(self):
        dag = DAG(inputs={}, steps=[], outputs={})
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf, 'workflow {')

    def test_outputspec_passthrough_still_emits_workflow_output_channel(self):
        dag = DAG(
            inputs={'x': Input('x', type='file?', desc=None, optional=True)},
            steps=[],
            outputs={'x': OutputSpec(type='file?', desc=None, optional=True, value=Input('x', type='file?', desc=None, optional=True))},
        )
        nf = transpile_dag_dict(dag.to_dict())
        self._assert_nf_contains(nf, 'x = x_ch')
        self.assertNotIn('emit:', nf)


    # Q5e: record-binding error messages -----------------------------------

    def test_record_binding_error_includes_fields(self):
        bad = {
            'inputs': {},
            'steps': [{
                'id': 'test',
                'type': 'task',
                'path': '/test.sh',
                'deps': [],
                'inputs': {'x': {'type': 'str'}},
                'bindings': {'x': {'source': 'record', 'fields': {'a': {'source': 'literal', 'value': 1}, 'b': {'source': 'literal', 'value': 2}}}},
                'outputs': {'out': {'type': 'file'}},
                'run': {},
                'script': 'echo hi',
            }],
            'outputs': {'out': {'type': 'file', 'desc': None, 'value': {'step': 'test', 'output': 'out'}}},
        }
        with self.assertRaisesRegex(ValueError, 'Record binding.*fields'):
            transpile_dag_dict(bad)

if __name__ == '__main__':
    unittest.main()
