#@  Stub merge task
#
# in
#   bcf [file] | input variant calls
#   outbase str | output base name
#
# out
#   bcf file = ${outbase}.bcf | output

cp "${bcf}" "${outbase}.bcf"
