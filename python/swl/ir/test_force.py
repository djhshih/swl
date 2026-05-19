import os
import unittest as ut

from swl.ir import node as ir
from swl.ir.dag import Literal, Merge, Record
from swl.ir.force import Forcer, ForceEnv, force_file
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
        }, root

    def test_force_saturated_workflow_produces_task_dag(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align', 'sort'])
        bam = data['steps'][1]['bindings']['bam']
        self.assertEqual(bam['source'], 'align')
        self.assertEqual(bam['output'], 'bam')
        self.assertIn('fastq1', data['inputs'])
        self.assertIn('outbase', data['inputs'])
        self.assertEqual(data['inputs']['fastq1']['type'], 'file')
        self.assertEqual(data['inputs']['fastq1']['desc'], 'fastq reads')
        self.assertEqual(data['inputs']['outbase']['type'], 'str')
        self.assertEqual(data['inputs']['outbase']['desc'], 'output base')
        self.assertIn('echo align', data['steps'][0]['script'])
        self.assertIn('outputs', data['steps'][0])
        self.assertEqual(data['steps'][0]['outputs']['bam']['default']['kind'], 'word')
        self.assertEqual(data['steps'][0]['deps'], [])
        self.assertEqual(data['steps'][1]['deps'], ['align'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bai', 'bam'])

    def test_force_partial_workflow_emits_task_with_remaining_inputs(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        data = dag.to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align'])
        self.assertEqual(sorted(data['inputs'].keys()), ['fastq1', 'fastq2', 'outbase'])
        self.assertEqual(sorted(data['outputs'].keys()), ['bam'])

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
        self.assertEqual(data['outputs']['bam']['step'], 'sort')
        self.assertEqual(data['outputs']['bam']['output'], 'bam')
        self.assertEqual(data['outputs']['bai']['step'], 'sort')
        self.assertEqual(data['outputs']['bai']['output'], 'bai')
        self.assertEqual(data['outputs']['bcf']['step'], 'call')
        self.assertEqual(data['outputs']['bcf']['output'], 'bcf')

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
        self.assertEqual(data['outputs']['bam']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['step'], 'align')

    def test_reused_computation_across_imported_workflow_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'shadow.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align', 'sub'])
        self.assertEqual(data['outputs']['bam']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['step'], 'sub')

    def test_partial_application_reuse_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'partial_reuse.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align'])
        self.assertEqual(data['outputs']['bam']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['step'], 'align')

    def test_workflow_partial_application_reuse_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'workflow_partial_reuse.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['mk'])
        self.assertEqual(data['outputs']['bam']['step'], 'mk')
        self.assertEqual(data['outputs']['bam2']['step'], 'mk')

    def test_nested_workflow_value_reuse_is_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'nested_lambda_reuse.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['sub'])
        self.assertEqual(data['outputs']['bam']['step'], 'sub')
        self.assertEqual(data['outputs']['bam2']['step'], 'sub')

    def test_partial_application_with_different_args_is_not_deduped(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'partial_different.swl'), files).to_dict()
        self.assertEqual([task['id'] for task in data['steps']], ['align', 'align_2'])
        self.assertEqual(data['outputs']['bam']['step'], 'align')
        self.assertEqual(data['outputs']['bam2']['step'], 'align_2')

    def test_map_force_produces_mapped_step_and_array_field(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'batch.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['align', 'merge'])
        self.assertIn('map', data['steps'][0])
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'input')
        self.assertEqual(data['steps'][1]['bindings']['bam']['source'], 'align')
        self.assertEqual(data['steps'][1]['deps'], ['align'])

    def test_map_lambda_forces_as_generated_mapped_workflow(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_lambda.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['map_lambda_1', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][1]['bindings']['bam']['source'], 'map_lambda_1')
        self.assertIn('x', data['inputs'])

    def test_map_lambda_with_inner_task_forces_as_generated_mapped_workflow(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_lambda_task.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['map_lambda_1', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        definition = data['steps'][0]['definition']
        self.assertEqual(definition['class'], 'Workflow')
        self.assertEqual([step['id'] for step in definition['dag']['steps']], ['sub'])
        self.assertEqual(data['steps'][1]['bindings']['bam']['source'], 'map_lambda_1')
        self.assertEqual(sorted(data['inputs'].keys()), ['x'])

    def test_map_partial_task_application_produces_mapped_step(self):
        files, root = self._files()
        lowered = lower_file(os.path.join(root, 'map_partial.swl'), files)
        mapped = lowered.body.bindings[1].value
        self.assertIsInstance(mapped.function, ir.Function)
        self.assertEqual(mapped.function.kind, 'workflow')
        data = force_file(os.path.join(root, 'map_partial.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['map_partial_1', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'input')
        self.assertEqual(data['steps'][1]['bindings']['bam']['source'], 'map_partial_1')

    def test_map_imported_workflow_produces_mapped_step(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_workflow.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['mk', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][0]['map']['ports'], ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'])
        self.assertEqual(sorted(data['inputs'].keys()), ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'])
        self.assertIn('fastq1', data['inputs'])

    def test_map_imported_workflow_with_partial_task_inside_produces_mapped_step(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'map_workflow_partial.swl'), files).to_dict()
        self.assertEqual([step['id'] for step in data['steps']], ['mkp', 'merge'])
        self.assertEqual(data['steps'][0]['type'], 'workflow')
        self.assertEqual(data['steps'][0]['map']['ports'], ['fastq1', 'fastq2', 'outbase'])

    def test_force_rejects_unnormalized_map_callable(self):
        with self.assertRaisesRegex(ValueError, 'map requires normalized executable callable during forcing'):
            Forcer().force_value(
                ir.Map(
                    ir.Lambda('x', ir.Block([], ir.Record({'bam': ir.Field(ir.Name('x'), 'bam')}))),
                    ir.Name('xs'),
                ),
                ForceEnv(),
            )

    def test_force_rejects_unsupported_ir_node(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported IR node during forcing: object'):
            Forcer().force_value(object(), None)

    def test_merge_of_records_collapses_during_force(self):
        value = Forcer().force_value(
            ir.Update(
                ir.Record({'a': ir.Literal(1)}),
                ir.Record({'b': ir.Literal(2)}),
            ),
            None,
        )
        self.assertEqual(value.fields['a'].value, 1)
        self.assertEqual(value.fields['b'].value, 2)

    def test_field_projection_prefers_right_side_of_merge(self):
        value = Forcer().force_value(
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

    def test_dependency_extraction_walks_merged_record_inputs(self):
        files, root = self._files()
        data = force_file(os.path.join(root, 'chain.swl'), files).to_dict()
        self.assertEqual(data['steps'][1]['deps'], ['align'])
        self.assertEqual(data['steps'][2]['deps'], ['sort'])

    def test_output_flattening_handles_nested_merged_records(self):
        outputs = Forcer()._final_outputs(
            Merge(
                Record({'a': Literal(1)}),
                Merge(Record({'b': Literal(2)}), Record({'c': Literal(3)})),
            )
        )
        self.assertEqual(sorted(outputs.keys()), ['a', 'b', 'c'])


if __name__ == '__main__':
    ut.main()
