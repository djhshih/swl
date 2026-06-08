#@  Stub sort task
#
# in
#   bam file | input bam
#   outbase str | output base name
#
# out
#   bam file = ${outbase}.bam | output
#   bai file = ${outbase}.bai | index

cp "${bam}" "${outbase}.bam"
cp "${bam}" "${outbase}.bai"
