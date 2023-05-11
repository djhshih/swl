## Align paired-end sequencing reads
#
# fun:
#   fastq1, fastq2, ref, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa -> bam
#
# in:
#   fastq1, fastq2 file | paired-end reads
#   ref file | reference sequence
#   ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa file
#     | reference bwa index files
#   outbase str | output base name
#
# out:
#   bam file = ${outbase}.bam
#     | output alignment
#
# run:
#   cpu int = 2
#   image size = djhshih/seqkit:0.1

bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam
