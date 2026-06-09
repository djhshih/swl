#@  Stub call task
#
# in
#   bam file | input bam
#   ref file | reference
#   ref_fai file | reference index
#   outbase str | output base name
#
# out
#   bcf file = ${outbase}.bcf | output

cat "${bam}" "${ref}" "${ref_fai}" > "${outbase}.bcf"
