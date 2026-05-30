import json
import os
import tempfile
import unittest as ut

from swl.ir.force import DAG, force_file
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
        restored = DAG.from_dict(data)
        self.assertEqual(restored.to_dict(), data)

    def test_partial_map_root_round_trip(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'map_root.swl'), files)
        data = dag.to_dict()
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


if __name__ == '__main__':
    ut.main()
