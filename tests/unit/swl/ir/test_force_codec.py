import json
import os
import tempfile
import unittest as ut

from swl.ir.dag import Field, Input, Literal, OutputSpec, StepCall
from swl.ir.force import DAG, Forcer, force_file
from swl.ir.lower import Lowerer


_ALIGN = '''# @ Align
# in
#   fastq1, fastq2 file
#   ref file
#   ref_fai file
#   outbase str
# out
#   bam file = ${outbase}.bam
# run
#   cpu = 2
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

_BAD_OUTBASE = '''# @ Bad
# out
#   outbase file = result.txt
echo bad
'''


class TestForceCodec(ut.TestCase):
    def _files(self):
        root = os.path.abspath('/virtual-force-codec')
        return {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sort.sh'): _SORT,
            os.path.join(root, 'bad_outbase.sh'): _BAD_OUTBASE,
            os.path.join(root, 'pipe.swl'): _PIPE,
            os.path.join(root, 'partial.swl'): '''align = import "align.sh"
align_hg38 = align {
  ref: "hg38.fa",
  ref_fai: "hg38.fa.fai",
  cpu: 4
}
align_hg38
''',
            os.path.join(root, 'bad_pipe.swl'): '''align = import "bad_outbase.sh"
sort = import "sort.sh"
align | sort
''',
            os.path.join(root, 'merge.sh'): '''# @ Merge
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
            os.path.join(root, 'map_root.swl'): '''call_variant = import "pipe.swl"
map call_variant
''',
            os.path.join(root, 'map_by_root.swl'): '''call_variant = import "pipe.swl"
map_by call_variant "outbase"
''',
        }, root

    def test_dag_round_trip_dict(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        restored = DAG.from_dict(data)
        self.assertEqual(restored.to_dict(), data)

    def test_dag_round_trip_json(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        payload = json.dumps(dag.to_dict())
        restored = DAG.from_dict(json.loads(payload))
        self.assertEqual(restored.to_dict(), dag.to_dict())

    def test_forced_dag_json_contains_no_legacy_array_field_encoding(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch.swl'), files)
        payload = json.dumps(dag.to_dict(), sort_keys=True)
        self.assertNotIn('mapped_value', payload)
        self.assertNotIn('array_field', payload)
        self.assertNotIn('"source": "step"', payload)

    def test_dag_write_and_read(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        td = tempfile.TemporaryDirectory()
        self.addCleanup(td.cleanup)
        path = os.path.join(td.name, 'plan.json')
        dag.write(path)
        restored = DAG.read(path)
        self.assertEqual(restored.to_dict(), dag.to_dict())

    def test_pipe_and_function_sources_compile_identically(self):
        root = os.path.abspath('/virtual-force-equivalence')
        call = '''# @ Call
# in
#   bam file
#   ref file
#   ref_fai file
#   outbase str
# out
#   bcf file = ${outbase}.bcf
#   | called variants
echo call
'''
        files = {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sort.sh'): _SORT,
            os.path.join(root, 'call.sh'): call,
            os.path.join(root, 'pipe.swl'): '''align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"

align | sort | call
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
        }
        pipe = force_file(os.path.join(root, 'pipe.swl'), files).to_dict()
        function = force_file(os.path.join(root, 'function.swl'), files).to_dict()
        self.assertEqual(function, pipe)

    def test_mapped_step_round_trip(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'batch.swl'), files)
        data = dag.to_dict()
        self.assertEqual(data['steps'][0]['map']['source']['source'], 'table')
        self.assertEqual(sorted(data['steps'][0]['map']['source']['columns'].keys()), ['fastq1', 'fastq2', 'outbase', 'ref', 'ref_fai'])
        restored = DAG.from_dict(data)
        self.assertEqual(restored.to_dict(), data)

    def test_partial_map_root_round_trip(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'map_root.swl'), files)
        data = dag.to_dict()
        self.assertEqual(data['steps'][0]['map']['source'], {'source': 'input', 'name': 'xs'})
        restored = DAG.from_dict(data)
        self.assertEqual(restored.to_dict(), data)

    def test_partial_map_by_root_round_trip_preserves_group_by(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'map_by_root.swl'), files)
        data = dag.to_dict()
        self.assertEqual(data['steps'][0]['map']['source'], {'source': 'input', 'name': 'xs'})
        self.assertEqual(data['steps'][0]['map']['group_by'], 'outbase')
        restored = DAG.from_dict(data)
        self.assertEqual(restored.to_dict(), data)

    def test_run_value_round_trip(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'partial.swl'), files)
        data = dag.to_dict()
        self.assertNotIn('default', data['steps'][0]['run']['cpu'])
        self.assertEqual(data['steps'][0]['run']['cpu']['value'], 4)
        restored = DAG.from_dict(data)
        self.assertEqual(restored.to_dict(), data)

    def test_outputs_round_trip_as_outputspec(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'pipe.swl'), files)
        data = dag.to_dict()
        self.assertIn('type', data['outputs']['bam'])
        self.assertIn('value', data['outputs']['bam'])
        self.assertEqual(data['outputs']['bam']['type'], 'file')
        self.assertEqual(data['outputs']['bam']['value'], {'source': 'step_output', 'step': 'sort', 'output': 'bam'})
        restored = DAG.from_dict(data)
        self.assertEqual(restored.to_dict(), data)

    def test_dag_validate_detects_circular_dependency(self):
        dag = DAG(
            inputs={},
            steps=[
                StepCall(id='a', path='/a.sh', bindings={}, outputs=[], task={'body': '', 'inputs': {}, 'outputs': {}, 'run': {}}, deps=['b']),
                StepCall(id='b', path='/b.sh', bindings={}, outputs=[], task={'body': '', 'inputs': {}, 'outputs': {}, 'run': {}}, deps=['a']),
            ],
            outputs={},
        )
        with self.assertRaisesRegex(ValueError, 'Circular dependency detected in DAG involving step: '):
            dag.validate()

    def test_dag_validate_detects_self_loop(self):
        dag = DAG(
            inputs={},
            steps=[
                StepCall(id='a', path='/a.sh', bindings={}, outputs=[], task={'body': '', 'inputs': {}, 'outputs': {}, 'run': {}}, deps=['a']),
            ],
            outputs={},
        )
        with self.assertRaisesRegex(ValueError, 'Circular dependency detected in DAG involving step: a'):
            dag.validate()

    def test_dag_validate_detects_unknown_dependency(self):
        dag = DAG(
            inputs={},
            steps=[
                StepCall(id='a', path='/a.sh', bindings={}, outputs=[], task={'body': '', 'inputs': {}, 'outputs': {}, 'run': {}}, deps=['nonexistent']),
            ],
            outputs={},
        )
        with self.assertRaisesRegex(ValueError, 'depends on unknown step'):
            dag.validate()

    def test_dag_from_dict_rejects_unknown_binding_source(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported binding source'):
            DAG.from_dict({
                'inputs': {},
                'steps': [],
                'outputs': {'x': {'source': 'mystery'}},
            })

    def test_lower_file_rejects_semantic_errors_before_force(self):
        files, root = self._files()
        with self.assertRaisesRegex(ValueError, 'Type mismatch for "outbase": file -> str'):
            Lowerer(files=files).lower_file(os.path.join(root, 'bad_pipe.swl'))

    # P4b: mapped-port validation --------------------------------------------

    def _mapped_step(self, id='m', scatter=None, broadcast=None, schema=None, extra_inputs=None):
        schema = schema or {'x': 'str'}
        input_names = list(schema.keys())
        m = StepCall(
            id=id, path=f'/{id}.sh', bindings={}, outputs=['r'],
            map={'source': {'source': 'input', 'name': 'xs'},
                 **({'scatter': scatter} if scatter is not None else {}),
                 **({'broadcast': broadcast} if broadcast is not None else {})},
            input_schema=schema, output_schema={'r': 'str'},
            task={'body': '', 'inputs': {n: {'type': schema[n]} for n in input_names},
                  'outputs': {'r': {'type': 'str'}}, 'run': {}},
        )
        return m

    def test_validate_mapped_step_missing_scatter(self):
        dag = DAG(
            inputs={'xs': Input('xs', '[str]')},
            steps=[self._mapped_step(scatter=None, broadcast=[])],
            outputs={},
        )
        with self.assertRaisesRegex(ValueError, 'missing map.scatter'):
            dag.validate()

    def test_validate_mapped_step_missing_broadcast(self):
        dag = DAG(
            inputs={'xs': Input('xs', '[str]')},
            steps=[self._mapped_step(scatter=['x'], broadcast=None)],
            outputs={},
        )
        with self.assertRaisesRegex(ValueError, 'missing map.broadcast'):
            dag.validate()

    def test_validate_mapped_step_input_in_both_ports_raises(self):
        dag = DAG(
            inputs={'xs': Input('xs', '[str]')},
            steps=[self._mapped_step(scatter=['x'], broadcast=['x'])],
            outputs={},
        )
        with self.assertRaisesRegex(ValueError, 'appears in both'):
            dag.validate()

    def test_validate_mapped_step_input_in_neither_port_raises(self):
        dag = DAG(
            inputs={'xs': Input('xs', '[str]')},
            steps=[self._mapped_step(scatter=[], broadcast=[], schema={'x': 'str', 'y': 'str'})],
            outputs={},
        )
        with self.assertRaisesRegex(ValueError, 'neither scatter nor broadcast'):
            dag.validate()

    def test_validate_mapped_step_valid_inputs_pass(self):
        dag = DAG(
            inputs={'xs': Input('xs', '[str]'), 'y': Input('y', 'str')},
            steps=[self._mapped_step(scatter=['x'], broadcast=['y'], schema={'x': 'str', 'y': 'str'})],
            outputs={},
        )
        dag.validate()  # should not raise


if __name__ == '__main__':
    ut.main()
