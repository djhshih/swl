#@ Concatenate multiple files
#
# in
#   files [file]   | files to concatenate
#   outname str
# out
#   combined file = ${outname}  | combined file

cat ${files} > ${combined}
