#@  Stub sort task
#
# in
#   bam file | input bam
#   outbase str | output base name
#
# out
#   bam file = ${outbase}.sorted.bam | output
#   bai file = ${outbase}.bai | index

cp "${bam}" "${outbase}.sorted.bam"
cp "${bam}" "${outbase}.bai"
