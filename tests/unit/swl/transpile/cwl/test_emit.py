import os
import unittest

from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.dag.forcer import force_file
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

    def test_transpile_function_workflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        self.assertEqual(cwl['cwlVersion'], 'v1.2')
        tools = [item for item in cwl['$graph'] if item['class'] == 'CommandLineTool']
        self.assertEqual(len(tools), 3)
        workflow = cwl['$graph'][-1]
        self.assertEqual(len(workflow['outputs']), 3)
        align = tools[0]
        self.assertEqual(align['baseCommand'], ['bash', 'script.sh'])

    def test_output_glob_uses_cwl_expression(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        align = next(item for item in cwl['$graph'] if item.get('id') == '#align')
        self.assertEqual(align['outputs'][0]['outputBinding']['glob'], "$(inputs.outbase + '.bam')")

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
            'steps': [{
                'id': 'align', 'path': '/tmp/align.sh', 'deps': [],
                'inputs': {'x': {'type': 'str', 'desc': None}},
                'bindings': {'x': {'source': 'merge', 'left': {'source': 'input', 'name': 'a'}, 'right': {'source': 'input', 'name': 'b'}}},
                'outputs': {'bam': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'x.bam'}]}, 'desc': None}},
                'run': {}, 'script': 'echo hi\n',
            }],
            'outputs': {'bam': {'type': 'file', 'desc': None, 'value': {'step': 'align', 'output': 'bam'}}},
        }
        with self.assertRaisesRegex(ValueError, 'contains a Merge value'):
            transpile_dag_dict(bad)

    def test_batch_mapped_task_emits_scatter_and_tab_column_input_type(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        align = next(step for step in workflow['steps'] if step['id'] == '#main/align')
        self.assertEqual(align['scatterMethod'], 'dotproduct')
        self.assertIn('scatter', align)
        self.assertGreater(len(align['scatter']), 0)
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
        self.assertIn('scatter', mk)
        self.assertGreater(len(mk['scatter']), 0)
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

    def test_map_by_transpile_emits_grouping_and_wrapper(self):
        dag_dict = {
            'inputs': {'sample': {'type': '[str]', 'desc': None}},
            'steps': [{
                'id': 'grouped',
                'type': 'workflow',
                'path': '/tmp/grouped.swl',
                'map': {'source': {'source': 'table', 'name': 'table', 'columns': {'sample': {'source': 'input', 'name': 'sample'}}}, 'group_by': 'sample', 'scatter': ['sample'], 'broadcast': []},
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
        cwl = transpile_dag_dict(dag_dict)
        self.assertEqual(cwl['cwlVersion'], 'v1.2')
        tool_ids = [item['id'] for item in cwl['$graph'] if item.get('class') in ('ExpressionTool', 'CommandLineTool')]
        self.assertIn('#group_grouped', tool_ids)
        self.assertIn('#wrap_grouped', tool_ids)
        workflow = cwl['$graph'][-1]
        step_ids = [s['id'] for s in workflow['steps']]
        self.assertIn('#main/group_grouped', step_ids)
        self.assertIn('#main/grouped', step_ids)
        wrap_step = next(s for s in workflow['steps'] if s['id'] == '#main/grouped')
        self.assertIn('scatter', wrap_step)
        self.assertEqual(wrap_step['scatterMethod'], 'dotproduct')

    def test_map_by_from_compiled_dag_emits_correct_structure(self):
        files, root = self._files()
        files[os.path.join(root, 'map_by_test.swl')] = (
            'tool = import "align.sh"\n'
            'map_by tool "outbase"\n'
        )
        dag = force_file(os.path.join(root, 'map_by_test.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        self.assertEqual(cwl['cwlVersion'], 'v1.2')
        group_tools = [item for item in cwl['$graph'] if item.get('class') == 'CommandLineTool' and 'group' in item.get('id', '')]
        self.assertTrue(group_tools, 'Should have at least one CommandLineTool for grouping')
        wrap_tools = [item for item in cwl['$graph'] if item.get('class') == 'CommandLineTool' and 'wrap' in item.get('id', '')]
        self.assertTrue(wrap_tools, 'Should have wrapper tools')
        workflow = cwl['$graph'][-1]
        self.assertTrue(any('scatter' in s for s in workflow['steps']), 'Should have scattered steps')

    def test_root_partial_map_transpiles_as_scattered_subworkflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'map_root.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        step = next(step for step in workflow['steps'] if step['id'] == '#main/call_variant')
        self.assertEqual(step['run'], '#call_variant')
        self.assertEqual(step['scatterMethod'], 'dotproduct')
        self.assertIn('scatter', step)
        self.assertGreater(len(step['scatter']), 0)
        outputs = {item['id']: item for item in workflow['outputs']}
        self.assertEqual(outputs['#main/bam']['type'], {'type': 'array', 'items': 'File'})
        subwf = next(item for item in cwl['$graph'] if item.get('id') == '#call_variant')
        self.assertEqual(subwf['class'], 'Workflow')

    def test_expr_interpolation_passthrough(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'bad_expr.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        bad_tool = next(item for item in cwl['$graph'] if item.get('id') == '#bad')
        bam_output = next(o for o in bad_tool['outputs'] if o['id'] == '#bad/bam')
        glob_expr = bam_output['outputBinding']['glob']
        self.assertIn('outbase', glob_expr)
        self.assertIn('/', glob_expr)
        req_classes = [r['class'] for r in bad_tool['requirements']]
        self.assertIn('InlineJavascriptRequirement', req_classes)

    def test_time_hint_emitted_for_task_with_time(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict())
        align = next(item for item in cwl['$graph'] if item.get('id') == '#align')
        self.assertIn('hints', align)
        time_hints = [h for h in align['hints'] if h.get('class') == 'TimeLimit']
        self.assertEqual(len(time_hints), 1)
        self.assertEqual(time_hints[0]['timeLimit'], 3150)

    def test_record_binding_step_input_emits_expression_tool(self):
        producer = StepCall(
            id='producer',
            path='/tmp/producer.sh',
            bindings={},
            outputs=['result'],
            task={
                'body': 'echo produce\n',
                'inputs': {},
                'outputs': {'result': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'out.txt'}]}, 'desc': None}},
                'run': {},
            },
        )
        consumer = StepCall(
            id='consumer',
            path='/tmp/consumer.sh',
            bindings={'input_x': Record({'a': Input('inp'), 'b': Field(producer, 'result')})},
            outputs=['out'],
            task={
                'body': 'echo consume\n',
                'inputs': {'input_x': {'type': 'str', 'desc': None}},
                'outputs': {'out': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'out2.txt'}]}, 'desc': None}},
                'run': {},
            },
            deps=['producer'],
        )
        dag = DAG(
            inputs={'inp': Input('inp', type='str', desc=None)},
            steps=[producer, consumer],
            outputs={'out': OutputSpec(type='file', desc=None, value=Field(consumer, 'out'))},
        )
        cwl = transpile_dag_dict(dag.to_dict())
        self.assertEqual(cwl['cwlVersion'], 'v1.2')
        expr_tools = [item for item in cwl['$graph'] if item.get('class') == 'ExpressionTool']
        self.assertEqual(len(expr_tools), 1, 'Should have exactly one ExpressionTool for record input binding')
        workflow = cwl['$graph'][-1]
        consumer_step = next(step for step in workflow['steps'] if step['id'] == '#main/consumer')
        rec_input = next(i for i in consumer_step['in'] if i['id'] == '#main/consumer/input_x')
        self.assertIn('rec_consumer_input_x', rec_input['source'])

    def test_optional_workflow_output_uses_null_union_type(self):
        dag = DAG(
            inputs={'x': Input('x', type='file?', desc=None, optional=True)},
            steps=[],
            outputs={'x': OutputSpec(type='file?', desc=None, optional=True, value=Input('x', type='file?', desc=None, optional=True))},
        )
        cwl = transpile_dag_dict(dag.to_dict())
        workflow = cwl['$graph'][-1]
        out = next(o for o in workflow['outputs'] if o['id'] == '#main/x')
        self.assertEqual(out['type'], ['null', 'File'])

    def test_nested_field_binding_uses_valueFrom(self):
        nest = StepCall(
            id='inner',
            path='/tmp/inner.sh',
            bindings={},
            outputs=['result'],
            task={
                'body': 'echo inner\n',
                'inputs': {},
                'outputs': {'result': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'out.txt'}]}, 'desc': None}},
                'run': {},
            },
        )
        leaf = StepCall(
            id='outer',
            path='/tmp/outer.sh',
            bindings={'x': Field(source=Field(source=nest, name='result'), name='value')},
            outputs=['out'],
            task={
                'body': 'echo outer\n',
                'inputs': {'x': {'type': 'str', 'desc': None}},
                'outputs': {'out': {'type': 'file', 'default': {'kind': 'word', 'parts': [{'kind': 'literal', 'text': 'out.txt'}]}, 'desc': None}},
                'run': {},
            },
            deps=['inner'],
        )
        dag = DAG(
            inputs={},
            steps=[nest, leaf],
            outputs={'out': OutputSpec(type='file', desc=None, value=Field(source=leaf, name='out'))},
        )
        cwl = transpile_dag_dict(dag.to_dict())
        outer_step = next(s for s in cwl['$graph'][-1]['steps'] if s['id'] == '#main/outer')
        x_input = next(i for i in outer_step['in'] if i['id'] == '#main/outer/x')
        self.assertEqual(x_input['source'], '#main/inner/result')
        self.assertEqual(x_input['valueFrom'], '$(value)')


    # ---- shell interpolation tests ----

    def test_shell_interpolation_resolves_inputs_and_run_params(self):
        from swl.transpile.cwl.emit import _interpolate_shell
        from swl.dag.node import StepCall

        step = StepCall(
            id='align', path='/t.sh', bindings={}, outputs=['bam'],
            task={
                'body': 'bwa mem -t ${cpu} ${ref} > ${outbase}.bam',
                'inputs': {'ref': {'type': 'file'}, 'outbase': {'type': 'str'}},
                'outputs': {'bam': {'type': 'file'}},
                'run': {'cpu': {'type': 'int', 'value': 2}},
            },
        )
        entry = _interpolate_shell('bwa mem -t ${cpu} ${ref} > ${outbase}.bam', step)
        self.assertIn('$(', entry)
        self.assertIn('"bwa mem -t "', entry)
        self.assertIn(' + 2 + ', entry)
        self.assertIn('inputs.ref.path', entry)
        self.assertIn('inputs.outbase', entry)

    def test_shell_interpolation_expression(self):
        from swl.transpile.cwl.emit import _interpolate_shell
        from swl.dag.node import StepCall

        step = StepCall(
            id='sort', path='/t.sh', bindings={}, outputs=['bam'],
            task={
                'body': 'samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam',
                'inputs': {'outbase': {'type': 'str'}},
                'outputs': {'bam': {'type': 'file'}},
                'run': {'cpu': {'type': 'int', 'value': 2}, 'memory': {'type': 'memory', 'value': 8192}},
            },
        )
        entry = _interpolate_shell(
            'samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam', step
        )
        self.assertIn('$(', entry)
        self.assertIn('(8192 / 2)', entry)
        self.assertIn('inputs.outbase', entry)

    def test_shell_interpolation_passes_through_builtins(self):
        from swl.transpile.cwl.emit import _interpolate_shell
        from swl.dag.node import StepCall

        step = StepCall(
            id='t', path='/t.sh', bindings={}, outputs=[],
            task={
                'body': 'echo $HOME',
                'inputs': {},
                'outputs': {},
                'run': {},
            },
        )
        entry = _interpolate_shell('echo $HOME', step)
        self.assertIn('"$HOME"', entry)
        self.assertIn('$(', entry)

    def test_shell_interpolation_skips_body_without_dollar(self):
        from swl.transpile.cwl.emit import _interpolate_shell
        from swl.dag.node import StepCall

        step = StepCall(
            id='t', path='/t.sh', bindings={}, outputs=[],
            task={
                'body': 'echo hello',
                'inputs': {},
                'outputs': {},
                'run': {},
            },
        )
        entry = _interpolate_shell('echo hello', step)
        self.assertEqual(entry, 'echo hello')

    def test_shell_interpolation_multi_line(self):
        from swl.transpile.cwl.emit import _interpolate_shell
        from swl.dag.node import StepCall

        step = StepCall(
            id='sort', path='/t.sh', bindings={}, outputs=['bam', 'bai'],
            task={
                'body': 'samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam\nsamtools index ${outbase}.bam ${outbase}.bai',
                'inputs': {'outbase': {'type': 'str'}},
                'outputs': {'bam': {'type': 'file'}, 'bai': {'type': 'file'}},
                'run': {'cpu': {'type': 'int', 'value': 2}, 'memory': {'type': 'memory', 'value': 8192}},
            },
        )
        entry = _interpolate_shell(
            'samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam\nsamtools index ${outbase}.bam ${outbase}.bai',
            step,
        )
        self.assertIn('\\n', entry)
        self.assertIn('$(', entry)
        self.assertIn('inputs.outbase', entry)

    def test_shell_interpolation_in_transpiled_tool(self):
        from swl.dag.node import DAG, Field, Input, OutputSpec, StepCall

        ref_input = Input('ref', type='file', desc=None)
        ref = Field(ref_input, 'ref')
        task = StepCall(
            id='align', path='/tmp/align.sh',
            bindings={'ref': ref, 'outbase': Field(Input('outbase'), 'outbase')},
            outputs=['bam'],
            task={
                'body': 'bwa mem -t ${cpu} ${ref} > ${outbase}.bam',
                'inputs': {
                    'ref': {'type': 'file', 'desc': None},
                    'outbase': {'type': 'str', 'desc': None},
                },
                'outputs': {
                    'bam': {'type': 'file', 'default': None, 'desc': None},
                },
                'run': {
                    'cpu': {'type': 'int', 'value': 2},
                },
            },
        )
        dag = DAG(
            inputs={
                'ref': Input('ref', type='file', desc=None),
                'outbase': Input('outbase', type='str', desc=None),
            },
            steps=[task],
            outputs={'bam': OutputSpec(type='file', desc=None, value=Field(task, 'bam'))},
        )
        cwl = transpile_dag_dict(dag.to_dict())
        tool = next(item for item in cwl['$graph'] if item.get('id') == '#align')
        entry = tool['requirements'][0]['listing'][0]['entry']
        self.assertIsInstance(entry, str)
        self.assertIn('$(', entry)
        self.assertIn('inputs.ref.path', entry)
        self.assertIn('inputs.outbase', entry)
        req_classes = [r['class'] for r in tool['requirements']]
        self.assertIn('InlineJavascriptRequirement', req_classes)

    def test_shell_interpolation_no_dollar_does_not_add_inlinejs(self):
        from swl.dag.node import DAG, Input, OutputSpec, StepCall, Field

        task = StepCall(
            id='t', path='/t.sh', bindings={}, outputs=['out'],
            task={
                'body': 'echo hello',
                'inputs': {},
                'outputs': {'out': {'type': 'file', 'default': None}},
                'run': {},
            },
        )
        dag = DAG(
            inputs={},
            steps=[task],
            outputs={'out': OutputSpec(type='file', desc=None, value=Field(task, 'out'))},
        )
        cwl = transpile_dag_dict(dag.to_dict())
        tool = next(item for item in cwl['$graph'] if item.get('id') == '#t')
        entry = tool['requirements'][0]['listing'][0]['entry']
        self.assertEqual(entry, 'echo hello')
        req_classes = [r['class'] for r in tool['requirements']]
        self.assertNotIn('InlineJavascriptRequirement', req_classes)


class TestCwlTypeHelpers(unittest.TestCase):
    def test_to_cwl_type_scalar(self):
        from swl.types import to_cwl_type
        self.assertEqual(to_cwl_type('str'), 'string')
        self.assertEqual(to_cwl_type('file'), 'File')
        self.assertEqual(to_cwl_type('int'), 'int')
        self.assertEqual(to_cwl_type('float'), 'float')

    def test_to_cwl_type_optional_scalar(self):
        from swl.types import to_cwl_type
        self.assertEqual(to_cwl_type('str?'), ['null', 'string'])
        self.assertEqual(to_cwl_type('file?'), ['null', 'File'])
        self.assertEqual(to_cwl_type('int?'), ['null', 'int'])
        self.assertEqual(to_cwl_type('float?'), ['null', 'float'])

    def test_to_cwl_type_array(self):
        from swl.types import to_cwl_type
        self.assertEqual(to_cwl_type('[str]'), {'type': 'array', 'items': 'string'})
        self.assertEqual(to_cwl_type('[file]'), {'type': 'array', 'items': 'File'})
        self.assertEqual(to_cwl_type('[int]'), {'type': 'array', 'items': 'int'})
        self.assertEqual(to_cwl_type('[float]'), {'type': 'array', 'items': 'float'})

    def test_to_cwl_type_optional_array(self):
        from swl.types import to_cwl_type
        self.assertEqual(to_cwl_type('[str]?'), ['null', {'type': 'array', 'items': 'string'}])
        self.assertEqual(to_cwl_type('[file]?'), ['null', {'type': 'array', 'items': 'File'}])

    def test_to_cwl_type_rejects_array_of_optional_items(self):
        from swl.types import to_cwl_type
        for invalid in ['[str?]', '[file?]', '[int?]', '[float?]']:
            with self.assertRaises(ValueError):
                to_cwl_type(invalid)

    def test_to_array_type_strips_optional_before_wrapping(self):
        from swl.types import to_array_type
        self.assertEqual(to_array_type('str'), '[str]')
        self.assertEqual(to_array_type('str?'), '[str]')
        self.assertEqual(to_array_type('file?'), '[file]')
        self.assertEqual(to_array_type('int?'), '[int]')
        self.assertEqual(to_array_type('float?'), '[float]')
        self.assertEqual(to_array_type('[str]'), '[str]')
        self.assertEqual(to_array_type('[str]?'), '[str]?')

    def test_normalize_swl_type_rejects_array_of_optional_items(self):
        from swl.types import normalize_swl_type
        for invalid in ['[str?]', '[file?]', '[int?]', '[float?]']:
            with self.assertRaises(ValueError):
                normalize_swl_type(invalid)


if __name__ == '__main__':
    unittest.main()
