import os
import unittest

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
        }, root

    def test_transpile_function_workflow(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict(), workflow_id='function')
        self.assertEqual(cwl['cwlVersion'], 'v1.0')
        self.assertEqual(cwl['$graph'][-1]['id'], '#function')
        tools = [item for item in cwl['$graph'] if item['class'] == 'CommandLineTool']
        self.assertEqual([item['id'] for item in tools], ['#align', '#sort', '#call'])
        workflow = cwl['$graph'][-1]
        outputs = {item['id']: item for item in workflow['outputs']}
        self.assertEqual(outputs['#function/bam']['outputSource'], '#function/t2/bam')
        self.assertEqual(outputs['#function/bai']['outputSource'], '#function/t2/bai')
        self.assertEqual(outputs['#function/bcf']['outputSource'], '#function/t3/bcf')
        align = tools[0]
        self.assertEqual(align['baseCommand'], ['bash', 'script.sh'])
        self.assertEqual(align['requirements'][0]['class'], 'InitialWorkDirRequirement')
        self.assertIn('bwa mem', align['requirements'][0]['listing'][0]['entry'])
        self.assertEqual(align['requirements'][1]['coresMin'], 2)
        self.assertEqual(align['requirements'][2]['dockerPull'], 'djhshih/seqkit:0.1')

    def test_output_glob_uses_cwl_expression(self):
        files, root = self._files()
        dag = force_file(os.path.join(root, 'function.swl'), files)
        cwl = transpile_dag_dict(dag.to_dict(), workflow_id='function')
        align = next(item for item in cwl['$graph'] if item.get('id') == '#align')
        self.assertEqual(align['outputs'][0]['outputBinding']['glob'], "$(inputs.outbase + '.bam')")


if __name__ == '__main__':
    unittest.main()
