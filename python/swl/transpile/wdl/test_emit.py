import os
import unittest

from swl.ir.dag import DAG, Field, Input, Literal, Merge, Record, StepCall
from swl.ir.force import force_file
from swl.transpile.wdl.emit import transpile_dag_dict


class TestWDLTranspile(unittest.TestCase):
    def _files(self):
        root = '/tmp/swl-wdl'
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
            os.path.join(root, 'map_root.swl'): '''call_variant = import "pipe.swl"

map call_variant
''',
        }, root

    def _assert_wdl_contains(self, wdl, *patterns):
        for pattern in patterns:
            self.assertIn(pattern, wdl, f'Expected WDL to contain: {pattern}')

    def test_transpile_function_workflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl,
            'version 1.1',
            'task align {',
            'task sort {',
            'task call {',
            'workflow main {',
            'input {',
        )
        self._assert_wdl_contains(wdl, 'File fastq1')
        self._assert_wdl_contains(wdl, 'File fastq2')
        self._assert_wdl_contains(wdl, 'File ref')
        self._assert_wdl_contains(wdl, 'String outbase')
        self._assert_wdl_contains(wdl, 'command <<<')
        self._assert_wdl_contains(wdl, 'bwa mem -t ~{cpu} ~{ref} ~{fastq1} ~{fastq2}')
        self._assert_wdl_contains(wdl, 'output {')
        self._assert_wdl_contains(wdl, 'File bam = sort.bam')
        self._assert_wdl_contains(wdl, 'File bai = sort.bai')
        self._assert_wdl_contains(wdl, 'File bcf = call.bcf')
        self._assert_wdl_contains(wdl, 'requirements {')
        self._assert_wdl_contains(wdl, 'cpu: 2')
        self._assert_wdl_contains(wdl, 'memory: "8192 MB"')
        self._assert_wdl_contains(wdl, 'container: "djhshih/seqkit:0.1"')

    def test_output_glob_uses_wdl_interpolation(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self.assertIn('~{outbase}.bam', wdl)

    def test_rejects_non_task_workflow_output(self):
        dag = DAG(inputs={}, steps=[], outputs={'x': Literal(1)})
        with self.assertRaisesRegex(ValueError, 'WDL does not support literal workflow outputs'):
            transpile_dag_dict(dag.to_dict())

    def test_partial_workflow_transpiles(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl, 'workflow main {')
        self._assert_wdl_contains(wdl, 'String outbase')
        self.assertIn('"hg38.fa"', wdl)

    def test_literal_top_level_workflow_is_rejected_before_transpile(self):
        files, root = self._files()
        with self.assertRaisesRegex(ValueError, 'Workflow must evaluate to a function'):
            force_file(os.path.join(root, 'literal_output.swl'), files)

    def test_imported_workflow_output_transpiles_as_workflow_step(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'import_partial.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl, 'workflow main_partial {')

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
            'outputs': {'bam': {'step': 'align', 'output': 'bam'}},
        }
        with self.assertRaisesRegex(ValueError, 'merge'):
            transpile_dag_dict(bad)


    def test_batch_mapped_task_emits_scatter(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl,
            'scatter (align_i in range(length(fastq1)))',
            'call align',
            'fastq1 = fastq1[align_i]',
            'Array[File] bam',
        )

    def test_batch_mapped_workflow_emits_scattered_subworkflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch_workflow.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl,
            'scatter',
            'Array[File] bam',
        )

    def test_batch_mapped_simple_lambda_emits_generated_scattered_subworkflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch_lambda.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl, 'scatter')

    def test_map_by_transpile_reports_explicit_grouping_not_supported(self):
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
                'definition': {'class': 'Workflow', 'dag': {'inputs': {'sample': {'type': 'str', 'desc': None}}, 'steps': [], 'outputs': {'sample': {'source': 'input', 'name': 'sample'}}}, 'inputs': {'sample': {'type': 'str', 'desc': None}}, 'outputs': {'sample': {'type': 'str'}}, 'body': '', 'run': {}},
            }],
            'outputs': {'sample': {'step': 'grouped', 'output': 'sample'}},
        }
        wdl = transpile_dag_dict(bad)
        self._assert_wdl_contains(wdl, 'collect_by_key')

    def test_root_partial_map_transpiles_as_scattered_subworkflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'map_root.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl, 'scatter')

    def test_expr_interpolation_in_output(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'bad_expr.swl'), files)
        wdl = transpile_dag_dict(dag.to_dict())
        self.assertIn('~{outbase / 2}', wdl)

    def test_outputspec_type_is_used_for_workflow_outputs(self):
        dag = DAG(
            inputs={'x': Input('x', type='file?', desc=None, optional=True)},
            steps=[],
            outputs={'x': type('OutputSpec', (), {'type': 'file?', 'desc': None, 'optional': True, 'value': Input('x', type='file?', desc=None, optional=True)})()},
        )
        wdl = transpile_dag_dict(dag.to_dict())
        self.assertIn('File? x = x', wdl)

    def test_wdl_type_mapping(self):
        from swl.transpile.wdl.emit import _wdl_type
        self.assertEqual(_wdl_type('file'), 'File')
        self.assertEqual(_wdl_type('str'), 'String')
        self.assertEqual(_wdl_type('int'), 'Int')
        self.assertEqual(_wdl_type('float'), 'Float')
        self.assertEqual(_wdl_type('[file]'), 'Array[File]')
        self.assertEqual(_wdl_type('[str]'), 'Array[String]')
        self.assertEqual(_wdl_type('[int]'), 'Array[Int]')
        self.assertEqual(_wdl_type('[float]'), 'Array[Float]')
        self.assertEqual(_wdl_type('file', optional=True), 'File?')
        self.assertEqual(_wdl_type('str', optional=True), 'String?')

    def test_task_name_sanitization(self):
        from swl.transpile.wdl.emit import _task_name
        self.assertEqual(_task_name('align'), 'align')
        self.assertEqual(_task_name('Align-BAM'), 'align_bam')
        self.assertEqual(_task_name('123align'), '_123align')
        self.assertEqual(_task_name('_align'), 'align')

    def test_bash_var_interpolation(self):
        from swl.transpile.wdl.emit import _interpolate_bash_vars
        self.assertEqual(
            _interpolate_bash_vars('bwa mem -t ${cpu} ${ref}'),
            'bwa mem -t ~{cpu} ~{ref}'
        )
        self.assertEqual(
            _interpolate_bash_vars('echo $HOME'),
            'echo $HOME'
        )
        self.assertEqual(
            _interpolate_bash_vars('no vars here'),
            'no vars here'
        )

    def test_empty_inputs_outputs(self):
        dag = DAG(
            inputs={},
            steps=[],
            outputs={},
        )
        wdl = transpile_dag_dict(dag.to_dict())
        self._assert_wdl_contains(wdl, 'version 1.1')
        self._assert_wdl_contains(wdl, 'workflow main {')

    def test_struct_emits_field_types_from_bindings(self):
        from swl.ir.dag import Input as DagInput, Literal, Record, StepCall, Field
        producer = StepCall(
            id='producer', path='/tmp/p.sh', bindings={}, outputs=['bam'],
            task={
                'body': 'echo', 'inputs': {},
                'outputs': {'bam': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'o.bam'}]}, 'desc': None}},
                'run': {},
            },
        )
        consumer = StepCall(
            id='consumer', path='/tmp/c.sh',
            bindings={
                'rec': Record({
                    'f_file': DagInput('inp_file'),
                    'f_str': DagInput('inp_str'),
                    'f_int': Literal(42),
                    'f_step': Field(producer, 'bam'),
                }),
            },
            outputs=['out'],
            task={
                'body': 'echo', 'inputs': {'rec': {'type': 'str', 'desc': None}},
                'outputs': {'out': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'o.txt'}]}, 'desc': None}},
                'run': {},
            },
            deps=['producer'],
        )
        dag = DAG(
            inputs={
                'inp_file': DagInput('inp_file', type='file', desc=None),
                'inp_str': DagInput('inp_str', type='str', desc=None),
            },
            steps=[producer, consumer],
            outputs={'out': Field(consumer, 'out')},
        )
        wdl = transpile_dag_dict(dag.to_dict())
        self.assertIn('struct _Rec_F_file_F_int_F_step_F_str {', wdl)
        self.assertIn('File f_file', wdl)
        self.assertIn('String f_str', wdl)
        self.assertIn('Int f_int', wdl)
        self.assertIn('File f_step', wdl)

    def test_none_literal_raises_error(self):
        from swl.transpile.wdl.emit import _literal_to_wdl
        with self.assertRaisesRegex(ValueError, 'None/null literals'):
            _literal_to_wdl(None)


if __name__ == '__main__':
    unittest.main()
