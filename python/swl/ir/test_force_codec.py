import json
import os
import tempfile
import unittest as ut

from swl.ir.force import DAG, force_file


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


class TestForceCodec(ut.TestCase):
    def _files(self):
        root = os.path.abspath('/virtual-force-codec')
        return {
            os.path.join(root, 'align.sh'): _ALIGN,
            os.path.join(root, 'sort.sh'): _SORT,
            os.path.join(root, 'pipe.swl'): _PIPE,
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


if __name__ == '__main__':
    ut.main()
