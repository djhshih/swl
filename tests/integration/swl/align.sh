#@  Stub align task
#
# in
#   fastq1 file | read 1
#   fastq2 file | read 2
#   ref file | reference
#   ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa file | ref index
#   outbase str | output base name
#
# out
#   bam file = ${outbase}.bam | output

cp "${fastq1}" "${outbase}.bam"
