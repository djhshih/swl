#@  Stub merge task
#
# in
#   bcf [file] | input variant calls
#   outbase str | output base name
#
# out
#   bcf file = ${outbase}.bcf | output

cat ${bcf} > "${outbase}.bcf"
