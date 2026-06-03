import os
import unittest as ut

import swl.ir.node as ir
from swl.dag.node import DAG, Field, Input, Literal, Merge, OutputSpec, Record, StepCall
from swl.dag.context import ForceEnv
from swl.dag.finalize import (
    _build_output_specs,
    _final_outputs,
    _finalize_dag,
    _flatten_merge_value,
    _flatten_outputs,
    _flatten_step_bindings,
    _infer_output_type,
    _prune_unused_inputs,
)
from swl.dag.forcer import force_file, make_force_state
from swl.dag.evaluator import force_value
from swl.ir.lower import lower_file


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
# run
#   cpu = 2
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
align_hg38 = align {
    ref: "hg38.fa",
    ref_fai: "hg38.fa.fai",
    cpu: 4
}
align_hg38
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

_FUNCTION = '''align = import "align.sh"
sort = import "sort.sh"
call = import "call.sh"

\\x ->
    a = align x
    s = sort ( x // a )
    c = call ( x // a // s )
    a // s // c
'''

_REUSE = '''align = import "align.sh"
\\x ->
    a = align x
    { bam: a.bam, bam2: a.bam }
'''

_SHADOW = '''align = import "align.sh"
sub = import "sub.swl"
\\x ->
    a = align x
    b = sub x
    { bam: a.bam, bam2: b.bam2 }
'''

_PARTIAL_REUSE = '''align = import "align.sh"
\\x ->
    f = align { ref: x.ref, ref_fai: x.ref_fai }
    a = f { fastq1: x.fastq1, fastq2: x.fastq2, outbase: x.outbase }
    b = f { fastq1: x.fastq1, fastq2: x.fastq2, outbase: x.outbase }
    { bam: a.bam, bam2: b.bam }
'''

_WORKFLOW_PARTIAL_REUSE = '''mk = import "mk_align.swl"
\\x ->
    f = mk { ref: x.ref, ref_fai: x.ref_fai }
    a = f { fastq1: x.fastq1, fastq2: x.fastq2, outbase: x.outbase }
    b = f { fastq1: x.fastq1, fastq2: x.fastq2, outbase: x.outbase }
    { bam: a.bam, bam2: b.bam }
'''

_NESTED_LAMBDA_REUSE = '''sub = import "sub.swl"
\\x ->
    outer = sub
    r1 = outer x
    r2 = outer x
    { bam: r1.bam2, bam2: r2.bam2 }
'''

_PARTIAL_DIFFERENT = '''align = import "align.sh"
\\x ->
    f = align { ref: x.ref, ref_fai: x.ref_fai }
    a = f { fastq1: x.fastq1, fastq2: x.fastq2, outbase: x.outbase }
    b = f { fastq1: x.fastq1, fastq2: x.fastq2, outbase: "other" }
    { bam: a.bam, bam2: b.bam }
'''

_MERGE = '''# @ Merge
# in
#   bam [file]
#   outbase str
# out
#   bam file = ${outbase}.bam
echo merge
'''

_BATCH = '''align = import "align.sh"
merge = import "merge.sh"
\\xs ->
    calls = map align xs
    merge { bam: calls.bam, outbase: "merged" }
'''

_MAP_LAMBDA = '''merge = import "merge.sh"
\\xs ->
    f = \\x -> { bam: x.bam }
    ys = map f xs
    merge { bam: ys.bam, outbase: "merged" }
'''

_MAP_LAMBDA_TASK = '''sub = import "sub.swl"
merge = import "merge.sh"
\\xs ->
    f = \\x -> sub x
    ys = map f xs
    merge { bam: ys.bam2, outbase: "merged" }
'''

_MAP_PARTIAL = '''align = import "align.sh"
merge = import "merge.sh"
\\xs ->
    f = align { ref: "hg38.fa", ref_fai: "hg38.fa.fai" }
    ys = map f xs
    merge { bam: ys.bam, outbase: "merged" }
'''

_MAP_WORKFLOW = '''mk = import "mk_align.swl"
merge = import "merge.sh"
\\xs ->
    ys = map mk xs
    merge { bam: ys.bam, outbase: "merged" }
'''

_MAP_WORKFLOW_PARTIAL = '''mkp = import "mk_align_partial.swl"
merge = import "merge.sh"
\\xs ->
    ys = map mkp xs
    merge { bam: ys.bam, outbase: "merged" }
'''

_MAP_ROOT = '''call_variant = import "pipe.swl"
map call_variant
'''

_MAP_BY_WORKFLOW = '''g = import "group_align.swl"
merge = import "merge.sh"
\\xs ->
    ys = map_by g "sample" xs
    merge { bam: ys.bam, outbase: "merged" }
'''


def _forcer(files=None):
    return make_force_state(files=files)


class TestForce(ut.TestCase):
    def _files(self):
        root = os.path.abspath('/virtual-force')
        return {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sort.sh'): _SORT,
            os.path.join(root, 'call.sh'): _CALL,
            os.path.join(root, 'pipe.swl'): _PIPE,
            os.path.join(root, 'partial.swl'): _PARTIAL,
            os.path.join(root, 'partial.swl'): _PARTIAL,
            os.path.join(root, 'chain.swl'): _CHAIN,
            os.path.join(root, 'function.swl'): _FUNCTION,
            os.path.join(root, 'reuse.swl'): _REUSE,
            os.path.join(root, 'sub.swl'): '''align = import "align.sh"
\\x ->
    a = align x
    { bam2: a.bam }
''',
            os.path.join(root, 'mk_align.swl'): '''align = import "align.sh"
\\x ->
    align x
''',
            os.path.join(root, 'group_align.swl'): '''align = import "align.sh"
\\x ->
    a = align x
    { sample: x.sample, bam: a.bam }
''',
            os.path.join(root, 'mk_align_partial.swl'): '''align = import "align.sh"
\\x ->
    f = align { ref: "hg38.fa", ref_fai: "hg38.fa.fai" }
    f x
''',
            os.path.join(root, 'shadow.swl'): _SHADOW,
            os.path.join(root, 'partial_reuse.swl'): _PARTIAL_REUSE,
            os.path.join(root, 'workflow_partial_reuse.swl'): _WORKFLOW_PARTIAL_REUSE,
            os.path.join(root, 'nested_lambda_reuse.swl'): _NESTED_LAMBDA_REUSE,
            os.path.join(root, 'partial_different.swl'): _PARTIAL_DIFFERENT,
            os.path.join(root, 'merge.sh'): _MERGE,
            os.path.join(root, 'batch.swl'): _BATCH,
            os.path.join(root, 'map_lambda.swl'): _MAP_LAMBDA,
            os.path.join(root, 'map_lambda_task.swl'): _MAP_LAMBDA_TASK,
            os.path.join(root, 'map_partial.swl'): _MAP_PARTIAL,
            os.path.join(root, 'map_workflow.swl'): _MAP_WORKFLOW,
            os.path.join(root, 'map_workflow_partial.swl'): _MAP_WORKFLOW_PARTIAL,
            os.path.join(root, 'map_root.swl'): _MAP_ROOT,
            os.path.join(root, 'map_by_workflow.swl'): _MAP_BY_WORKFLOW,
        }, root

    def test_force_saturated_workflow_produces_task_dag(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align', 'sort'])
        bam = data['steps'][1]['bindings']['bam']
        self.assertEqual(bam['step'], 'align')
        self.assertEqual(bam['output'], 'bam')
        self.assertEqual(data['steps'][1]['deps'], ['align'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bai', 'bam'])

    def test_force_partial_workflow_emits_task_with_remaining_inputs(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align'])
        self.assertEqual(sorted(data['inputs'].keys()), ['fastq1', 'fastq2', 'outbase'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bam'])

    def test_force_partial_map_root_materializes_batch_workflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'map_root.swl'), files)
        data = dag.to_dict()
        self.assertEqual(sorted(data['inputs'].keys()), ['xs'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bai', 'bam'])
        self.assertEqual([step['id'] for step in data['steps']], ['call_variant'])
        self.assertEqual(data['steps'][0]['map']['source'], {'source': 'input', 'name': 'xs'})
        self.assertEqual(sorted(data['steps'][0]['outputs'].keys()), ['bai', 'bam'])

    def test_task_run_value_prefers_partial_application_over_task_default(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'partial.swl'), files).to_dict()
        self.assertEqual(data['steps'][0]['run']['cpu']['type'], 'int')
        self.assertNotIn('default', data['steps'][0]['run']['cpu'])
        self.assertEqual(data['steps'][0]['run']['cpu']['value'], 4)

    def test_chain_root_is_instantiated_during_force(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'chain.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align', 'sort', 'call'])
        self.assertIn('fastq1', data['inputs'])
        self.assertIn('ref_fai', data['inputs'])
        self.assertEqual(data['inputs']['ref_fai']['type'], 'file')
        self.assertEqual(data['steps'][2]['deps'], ['sort'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bai', 'bam', 'bcf'])
        self.assertEqual(data['outputs']['bam']['type'], 'file')
        self.assertEqual(data['outputs']['bam']['value']['step'], 'sort')
        self.assertEqual(data['outputs']['bam']['value']['output'], 'bam')
        self.assertEqual(data['outputs']['bai']['value']['step'], 'sort')
        self.assertEqual(data['outputs']['bai']['value']['output'], 'bai')
        self.assertEqual(data['outputs']['bcf']['value']['step'], 'call')
        self.assertEqual(data['outputs']['bcf']['value']['output'], 'bcf')

    def test_serialized_dag_is_self_contained(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        text = repr(data)
        self.assertNotIn('Import(', text)
        self.assertNotIn('workflow', text)
        self.assertEqual(data['steps'][0]['path'], os.path.join(root, 'align.sh'))
        self.assertIn('script', data['steps'][0])
        self.assertIn('outputs', data['steps'][0])

    def test_function_and_chain_compile_to_same_shape(self):
        files, root = self._files()
        chain = force_file(os.path.join(root, 'chain.swl'), files).to_dict()
        function = force_file(os.path.join(root, 'function.swl'), files).to_dict()
        self.assertEqual(chain['inputs'], function['inputs'])
        self.assertEqual(chain['outputs'], function['outputs'])
        self.assertEqual(
            [(task['id'], task['deps'], sorted(task['bindings'].keys()), sorted(task['outputs'].keys())) for task in chain['steps']],
            [(task['id'], task['deps'], sorted(task['bindings'].keys()), sorted(task['outputs'].keys())) for task in function['steps']],
        )

    def test_reused_variable_forces_once(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'reuse.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align'])
        self.assertEqual(data['outputs']['bam']['value']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['value']['step'], 'align')

    def test_reused_computation_across_imported_workflow_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'shadow.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align', 'sub'])
        self.assertEqual(data['outputs']['bam']['value']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['value']['step'], 'sub')

    def test_partial_application_reuse_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'partial_reuse.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align'])
        self.assertEqual(data['outputs']['bam']['value']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['value']['step'], 'align')

    def test_workflow_partial_application_reuse_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'workflow_partial_reuse.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['mk'])
        self.assertEqual(data['outputs']['bam']['value']['step'], 'mk')
        self.assertEqual(data['outputs']['bam2']['value']['step'], 'mk')

    def test_nested_workflow_value_reuse_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'nested_lambda_reuse.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['sub'])
        self.assertEqual(data['outputs']['bam']['value']['step'], 'sub')
        self.assertEqual(data['outputs']['bam2']['value']['step'], 'sub')

    def test_partial_application_with_different_args_is_not_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'partial_different.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align', 'align_2'])
        self.assertEqual(data['outputs']['bam']['value']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['value']['step'], 'align_2')

    def test_map_force_produces_mapped_step_and_tab_column_binding(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'batch.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['align', 'merge'])
        self.assertIn('map', data['steps'][0])
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'table')
        self.assertEqual(sorted(data['steps'][0]['map']['source']['columns'].keys()), ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'])
        self.assertEqual(data['steps'][1]['bindings']['bam']['step'], 'align')
        self.assertEqual(data['steps'][1]['deps'], ['align'])

    def test_map_lambda_forces_as_generated_mapped_workflow(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_lambda.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['map_lambda_1', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][1]['bindings']['bam']['step'], 'map_lambda_1')
        self.assertIn('bam', data['inputs'])

    def test_map_lambda_with_inner_task_forces_as_generated_mapped_workflow(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_lambda_task.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['map_lambda_1', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        definition = data['steps'][0]['definition']
        self.assertEqual(definition['class'], 'Workflow')
        self.assertEqual([step['id'] for step in definition['dag']['steps']], ['sub'])
        self.assertEqual(data['steps'][1]['bindings']['bam']['step'], 'map_lambda_1')
        self.assertEqual(sorted(data['inputs'].keys()), ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'])

    def test_map_partial_task_application_produces_mapped_step(self):
        files, root = self._files()
        lowered = lower_file(os.path.join(root, 'map_partial.swl'), files)
        mapped = lowered.body.bindings[1].value
        self.assertIsInstance(mapped.function, ir.Function)
        self.assertEqual(mapped.function.kind, 'workflow')
        data = force_file(os.path.join(root, 'map_partial.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['map_partial_1', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'table')
        self.assertEqual(data['steps'][1]['bindings']['bam']['step'], 'map_partial_1')

    def test_map_imported_workflow_produces_mapped_step(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_workflow.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['mk', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'table')
        self.assertEqual(data['steps'][0]['input_schema'], {'fastq1': 'file', 'fastq2': 'file', 'outbase': 'str', 'ref': 'file', 'ref_fai': 'file'})
        self.assertEqual(sorted(data['inputs'].keys()), ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'])
        self.assertIn('fastq1', data['inputs'])

    def test_map_imported_workflow_with_partial_task_inside_produces_mapped_step(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_workflow_partial.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['mkp', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'table')
        self.assertEqual(sorted(data['steps'][0]['map']['source']['columns'].keys()), ['fastq1', 'fastq2', 'outbase'])
        self.assertEqual(data['steps'][0]['input_schema'], {'fastq1': 'file', 'fastq2': 'file', 'outbase': 'str'})

    def test_mapped_table_source_uses_explicit_logical_table_metadata(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'batch.swl'), files).to_dict()
        source = data['steps'][0]['map']['source']
        self.assertEqual(source['source'], 'table')
        self.assertEqual(source['name'], 'table')
        self.assertEqual(
            source['columns'],
            {
                'fastq1': {'source': 'input', 'name': 'fastq1'},
                'fastq2': {'source': 'input', 'name': 'fastq2'},
                'outbase': {'source': 'input', 'name': 'outbase'},
                'ref': {'source': 'input', 'name': 'ref'},
                'ref_fai': {'source': 'input', 'name': 'ref_fai'},
            },
        )

    def test_mapped_step_contains_scatter_and_broadcast_ports(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'batch.swl'), files).to_dict()
        map_info = data['steps'][0]['map']
        self.assertIn('scatter', map_info)
        self.assertIn('broadcast', map_info)
        self.assertEqual(
            sorted(map_info['scatter']),
            ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'],
        )
        self.assertEqual(map_info['broadcast'], [])

    def test_map_by_imported_workflow_preserves_grouping_metadata(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_by_workflow.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['g', 'merge'])
        self.assertEqual(data['steps'][0]['map']['group_by'], 'sample')
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'table')

    def test_force_rejects_unnormalized_map_callable(self):
        env = ForceEnv()
        env.bind('xs', Literal(42))
        with self.assertRaisesRegex(ValueError, 'map requires normalized executable callable during forcing'):
            force_value(
                _forcer(),
                ir.Map(
                    ir.Lambda('x', ir.Block([], ir.Record({'bam': ir.Field(ir.Name('x'), 'bam')}))),
                    ir.Name('xs'),
                ),
                env,
            )

    def test_force_rejects_unsupported_ir_node(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported IR node during forcing: object'):
            force_value(_forcer(), object(), None)

    def test_merge_of_records_collapses_during_force(self):
        value = force_value(
            _forcer(),
            ir.Update(
                ir.Record({'a': ir.Literal(1)}),
                ir.Record({'b': ir.Literal(2)}),
            ),
            None,
        )
        self.assertEqual(value.fields['a'].value, 1)
        self.assertEqual(value.fields['b'].value, 2)

    def test_field_projection_prefers_right_side_of_merge(self):
        value = force_value(
            _forcer(),
            ir.Field(
                ir.Update(
                    ir.Record({'a': ir.Literal(1)}),
                    ir.Record({'a': ir.Literal(2)}),
                ),
                'a',
            ),
            None,
        )
        self.assertEqual(value.value, 2)

    def test_table_update_on_step_call_raises_error(self):
        with self.assertRaisesRegex(ValueError, 'Record update.*task/workflow call result'):
            force_value(
                _forcer(),
                ir.Update(ir.Name('align'), ir.Record({'ref': ir.Literal('hg38.fa')})),
                ForceEnv(values={'align': StepCall('align', '/align.sh', {}, ['out'])}),
            )

    def test_table_update_right_side_step_call_raises_error(self):
        with self.assertRaisesRegex(ValueError, 'Record update.*task/workflow call result'):
            force_value(
                _forcer(),
                ir.Update(ir.Record({'ref': ir.Literal('hg38.fa')}), ir.Name('align')),
                ForceEnv(values={'align': StepCall('align', '/align.sh', {}, ['out'])}),
            )

    def test_table_update_input_record_raises_error(self):
        forcer = _forcer()
        forcer.inputs['xs'] = Input('xs', '[str]')
        with self.assertRaisesRegex(ValueError, 'Table-record update.*not supported'):
            force_value(
                forcer,
                ir.Update(ir.Input('xs'), ir.Record({'ref': ir.Literal('hg38.fa')})),
                ForceEnv(),
            )

    def test_table_update_record_input_raises_error(self):
        forcer = _forcer()
        forcer.inputs['xs'] = Input('xs', '[str]')
        with self.assertRaisesRegex(ValueError, 'Record-table update.*not supported'):
            force_value(
                forcer,
                ir.Update(ir.Record({'ref': ir.Literal('hg38.fa')}), ir.Input('xs')),
                ForceEnv(),
            )

    def test_table_update_right_side_step_call_raises_error_via_name(self):
        env = ForceEnv()
        step = StepCall(id='align', path='/a.sh', bindings={}, outputs=['bam'])
        env.bind('ys', step)
        with self.assertRaisesRegex(ValueError, 'Record update.*on a task/workflow call result'):
            force_value(
                _forcer(),
                ir.Update(
                    ir.Record({'extra': ir.Literal('x')}),
                    ir.Name('ys'),
                ),
                env,
            )

    def test_dependency_extraction_walks_merged_record_inputs(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'chain.swl'), files).to_dict()
        self.assertEqual(data['steps'][1]['deps'], ['align'])
        self.assertEqual(data['steps'][2]['deps'], ['sort'])

    def test_output_flattening_handles_nested_merged_records(self):
        outputs = _final_outputs(
            _forcer(),
            Merge(
                Record({'a': Literal(1)}),
                Merge(Record({'b': Literal(2)}), Record({'c': Literal(3)})),
            )
        )
        self.assertEqual(sorted(outputs.keys()), ['a', 'b', 'c'])

    def test_merge_flattening_collapses_record_merges_in_outputs(self):
        forcer = _forcer()
        forcer.inputs['x'] = Input('x')
        forcer.steps = []
        outputs = _flatten_outputs(forcer, {
            'result': Merge(
                Record({'a': Literal(1), 'b': Literal(2)}),
                Record({'c': Literal(3)}),
            )
        })
        self.assertEqual(sorted(outputs.keys()), ['a', 'b', 'c'])
        self.assertIsInstance(outputs['a'], Literal)

    def test_merge_flattening_collapses_input_record_merge_in_outputs(self):
        forcer = _forcer()
        forcer.inputs['x'] = Input('x')
        forcer.steps = []
        outputs = _flatten_outputs(forcer, {
            'result': Merge(
                Input('x'),
                Record({'y': Literal(1)}),
            )
        })
        self.assertEqual(sorted(outputs.keys()), ['y'])
        self.assertEqual(outputs['y'].value, 1)

    def test_merge_flattening_step_bindings(self):
        forcer = _forcer()
        step = StepCall(id='test', path='/test.sh', bindings={
            'input': Merge(
                Record({'x': Literal(1)}),
                Record({'y': Literal(2)}),
            ),
        }, outputs=['out'])
        forcer.steps = [step]
        _flatten_step_bindings(forcer)
        self.assertEqual(sorted(step.bindings.keys()), ['x', 'y'])

    def test_merge_flattening_rejects_incompatible_non_record(self):
        forcer = _forcer()
        with self.assertRaisesRegex(ValueError, 'Cannot flatten merge'):
            _flatten_merge_value(forcer,
                Merge(Literal(1), Literal(2)),
            )

    def test_merge_flattening_preserves_single_non_record(self):
        forcer = _forcer()
        result = _flatten_merge_value(forcer,
            Merge(
                Record({'a': Literal(1)}),
                Input('x'),
            )
        )
        self.assertIsInstance(result, Record)
        self.assertIn('a', result.fields)
        self.assertEqual(result.fields['a'].value, 1)

    # P4d: deeply nested merge trees -----------------------------------------

    def test_deeply_nested_merge_tree_flattens_to_fields(self):
        forcer = _forcer()
        result = _flatten_merge_value(forcer,
            Merge(
                Merge(
                    Record({'a': Literal(1)}),
                    Record({'b': Literal(2)}),
                ),
                Merge(
                    Record({'c': Literal(3)}),
                    Record({'d': Literal(4)}),
                ),
            )
        )
        self.assertIsInstance(result, Record)
        self.assertEqual(set(result.fields.keys()), {'a', 'b', 'c', 'd'})

    def test_deeply_nested_merge_in_step_bindings_vanishes_from_final_dag(self):
        forcer = _forcer()
        forcer.inputs['x'] = Input('x')
        step = StepCall(id='test', path='/test.sh', bindings={
            'in': Merge(
                Merge(
                    Record({'a': Literal(1)}),
                    Record({'b': Literal(2)}),
                ),
                Input('x'),
            ),
        }, outputs=['out'])
        forcer.steps = [step]
        _flatten_step_bindings(forcer)
        self.assertEqual(set(step.bindings.keys()), {'a', 'b'})

    def test_deeply_nested_merge_in_outputs_flattens_completely(self):
        forcer = _forcer()
        forcer.inputs['x'] = Input('x')
        forcer.steps = []
        flat = _flatten_outputs(forcer, {
            'result': Merge(
                Merge(
                    Record({'a': Literal(1)}),
                    Record({'b': Literal(2)}),
                ),
                Merge(
                    Record({'c': Literal(3)}),
                    Input('x'),
                ),
            ),
        })
        self.assertEqual(sorted(flat.keys()), ['a', 'b', 'c'])
        for name in ('a', 'b', 'c'):
            self.assertNotIsInstance(flat[name], (Merge, Record),
                                     f'{name} should be a leaf binding, not Merge/Record')

    def test_merge_flattening_three_way_records_merges_all_fields(self):
        forcer = _forcer()
        result = _flatten_merge_value(forcer,
            Merge(
                Merge(
                    Record({'a': Literal(1), 'b': Literal(2)}),
                    Record({'c': Literal(3)}),
                ),
                Record({'d': Literal(4), 'e': Literal(5)}),
            )
        )
        self.assertIsInstance(result, Record)
        self.assertEqual(set(result.fields.keys()), {'a', 'b', 'c', 'd', 'e'})

    def test_merge_flattening_overlapping_fields_right_wins(self):
        forcer = _forcer()
        result = _flatten_merge_value(forcer,
            Merge(
                Merge(
                    Record({'x': Literal(1)}),
                    Record({'x': Literal(2)}),
                ),
                Record({'x': Literal(3)}),
            )
        )
        self.assertEqual(result.fields['x'].value, 3)

    # P4c: optionality propagation -------------------------------------------

    def test_force_ir_input_optionality_propagates(self):
        forcer = _forcer()
        val = force_value(forcer, ir.Input('x', type='file?', desc=None), None)
        self.assertTrue(forcer.inputs['x'].optional)
        self.assertEqual(forcer.inputs['x'].type, 'file?')

    def test_force_ir_input_required_when_no_question_mark(self):
        forcer = _forcer()
        val = force_value(forcer, ir.Input('x', type='file', desc=None), None)
        self.assertFalse(forcer.inputs['x'].optional)

    def test_force_ir_input_required_when_type_none(self):
        forcer = _forcer()
        val = force_value(forcer, ir.Input('x', type=None, desc=None), None)
        self.assertFalse(forcer.inputs['x'].optional)

    def test_optionality_serializes_in_dag_when_true(self):
        dag = DAG(
            inputs={'x': Input('x', 'file?', desc=None, optional=True)},
            steps=[],
            outputs={'y': OutputSpec(type='file', value=Input('x', 'file?', desc=None, optional=True))},
        )
        d = dag.to_dict()
        self.assertEqual(d['inputs']['x']['optional'], True)

    def test_optionality_omitted_from_dag_when_false(self):
        dag = DAG(
            inputs={'x': Input('x', 'file', desc=None, optional=False)},
            steps=[],
            outputs={'y': OutputSpec(type='file', value=Literal(42))},
        )
        d = dag.to_dict()
        self.assertNotIn('optional', d['inputs']['x'])

    # P4a: non-saturating records in DAG ------------------------------------

    def test_non_saturating_record_in_step_binding_gets_flattened(self):
        forcer = _forcer()
        step = StepCall(id='test', path='/test.sh', bindings={
            'params': Record({'x': Literal(1), 'y': Literal(2)}),
        }, outputs=['out'],
            task={'body': '', 'inputs': {'x': {'type': 'int'}, 'y': {'type': 'int'}},
                  'outputs': {'out': {'type': 'int'}}, 'run': {}},
        )
        forcer.steps = [step]
        _flatten_step_bindings(forcer)
        self.assertNotIn('params', step.bindings)
        self.assertIn('x', step.bindings)
        self.assertIn('y', step.bindings)
        self.assertEqual(step.bindings['x'].value, 1)
        self.assertEqual(step.bindings['y'].value, 2)

    def test_non_saturating_record_in_outputs_gets_flattened(self):
        forcer = _forcer()
        forcer.steps = []
        flat = _flatten_outputs(forcer, {
            'nested': Record({
                'a': Literal(1),
                'b': Record({'c': Literal(3)}),
            }),
        })
        self.assertIn('a', flat)
        self.assertIn('b', flat)
        self.assertEqual(flat['a'].value, 1)

    # P1c: merge output patterns -------------------------------------------

    def test_merge_record_field_stepcall_preserves_field(self):
        step = StepCall(id='producer', path='/p.sh', bindings={}, outputs=['result'])
        forcer = _forcer()
        forcer.steps = []
        result = _flatten_merge_value(forcer,
            Merge(
                Record({'a': Literal(1)}),
                Field(step, 'result'),
            )
        )
        self.assertIsInstance(result, Record)
        self.assertIn('a', result.fields)
        self.assertIn('result', result.fields)
        self.assertIs(result.fields['result'].source, step)

    def test_merge_record_field_input_preserves_field(self):
        forcer = _forcer()
        forcer.inputs['x'] = Input('x', 'file')
        result = _flatten_merge_value(forcer,
            Merge(
                Record({'a': Literal(1)}),
                Field(Input('x'), 'field_x'),
            )
        )
        self.assertIsInstance(result, Record)
        self.assertIn('a', result.fields)
        self.assertIn('field_x', result.fields)

    def test_merge_record_with_bare_input_drops_input(self):
        forcer = _forcer()
        forcer.inputs['x'] = Input('x', 'file')
        result = _flatten_merge_value(forcer,
            Merge(
                Record({'a': Literal(1)}),
                Input('x'),
            )
        )
        self.assertIsInstance(result, Record)
        self.assertIn('a', result.fields)
        self.assertNotIn('x', result.fields)

    def test_end_to_end_merge_at_output_level_leaves_no_merge_in_dag(self):
        step = StepCall(
            id='producer', path='/p.sh', bindings={}, outputs=['out'],
            task={'body': '', 'inputs': {}, 'outputs': {'out': {'type': 'file'}}, 'run': {}},
        )
        forcer = _forcer()
        forcer.inputs['x'] = Input('x', 'file')
        forcer.steps = [step]
        dag = DAG(
            inputs=dict(forcer.inputs),
            steps=[step],
            outputs={},
        )
        outputs = _flatten_outputs(forcer, {
            'result': Merge(
                Record({'a': Literal(1)}),
                Record({'b': Field(step, 'out')}),
            ),
        })
        for name, value in outputs.items():
            self.assertNotIsInstance(
                value, (Merge, Record),
                f'Output {name!r} should not be Merge/Record after flattening, got {type(value).__name__}',
            )


    # Q5a: output type inference for complex bindings -----------------------

    def test_infer_output_type_nested_field_chain_resolves_terminal_step(self):
        step = StepCall(
            id='producer', path='/p.sh', bindings={}, outputs=['inner'],
            task={
                'body': '', 'inputs': {}, 'run': {},
                'outputs': {'inner': {'type': 'file', 'desc': 'inner file'}},
            },
        )
        forcer = _forcer()
        forcer.steps = [step]
        outer = Field(Field(step, 'inner'), 'leaf')
        typ = _infer_output_type(forcer, outer)
        self.assertEqual(typ, 'file')

    def test_infer_output_type_nested_field_chain_mapped_resolves_to_array(self):
        step = StepCall(
            id='producer', path='/p.sh', bindings={}, outputs=['inner'],
            task={
                'body': '', 'inputs': {}, 'run': {},
                'outputs': {'inner': {'type': 'file', 'desc': ''}},
            },
            map={'source': {'source': 'input', 'name': 'xs'}},
        )
        forcer = _forcer()
        forcer.steps = [step]
        outer = Field(Field(step, 'inner'), 'leaf')
        typ = _infer_output_type(forcer, outer)
        self.assertEqual(typ, '[file]')

    def test_infer_output_type_record_returns_none(self):
        forcer = _forcer()
        typ = _infer_output_type(forcer, Record({'a': Literal(1)}))
        self.assertIsNone(typ)

    # Q5b: output descriptions ---------------------------------------------

    def test_output_desc_propagated_from_step_output(self):
        step = StepCall(
            id='producer', path='/p.sh', bindings={}, outputs=['bam'],
            task={
                'body': '', 'inputs': {}, 'run': {},
                'outputs': {'bam': {'type': 'file', 'desc': 'aligned BAM'}},
            },
        )
        forcer = _forcer()
        forcer.inputs['x'] = Input('x', 'file')
        forcer.steps = [step]
        specs = _build_output_specs(forcer, {
            'result': Field(step, 'bam'),
        })
        self.assertEqual(specs['result'].desc, 'aligned BAM')

    def test_output_desc_none_for_input_source(self):
        forcer = _forcer()
        forcer.inputs['x'] = Input('x', 'file')
        specs = _build_output_specs(forcer, {
            'result': Input('x'),
        })
        self.assertIsNone(specs['result'].desc)

    def test_output_desc_none_for_literal_source(self):
        forcer = _forcer()
        specs = _build_output_specs(forcer, {
            'result': Literal(42),
        })
        self.assertIsNone(specs['result'].desc)

    # Q5c: optionality in prune_unused_inputs -----------------------------

    def test_map_by_source_input_created_through_input_method(self):
        step = StepCall(
            id='grouped', path='/g.sh', bindings={}, outputs=['out'],
            task={'body': '', 'inputs': {}, 'outputs': {'out': {'type': 'file'}}, 'run': {}},
            map={'source': {'source': 'input', 'name': 'xs'}, 'group_by': 'key'},
        )
        forcer = _forcer()
        forcer.steps = [step]
        forcer.inputs['xs'] = Input('xs', '[str]', optional=True)
        unused = set()
        _prune_unused_inputs(forcer, None, {})
        self.assertIn('xs', forcer.inputs)

    # Q5e: compile-time guard for Record in outputs -----------------------

    def test_flatten_outputs_handles_nested_record_via_final_outputs(self):
        step = StepCall(
            id='producer', path='/p.sh', bindings={}, outputs=['out'],
            task={'body': '', 'inputs': {}, 'outputs': {'out': {'type': 'file'}}, 'run': {}},
        )
        forcer = _forcer()
        forcer.inputs['x'] = Input('x', 'file')
        forcer.steps = [step]
        dag = _finalize_dag(forcer,
            Record({'nested': Record({'a': Literal(1)})})
        )
        self.assertIn('a', dag.outputs)
        self.assertNotIn('nested', dag.outputs)


if __name__ == '__main__':
    ut.main()
