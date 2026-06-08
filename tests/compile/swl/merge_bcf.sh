#@  Merge mutations
# in
#   bcf      [file]                  | input mutation calls
#   outbase  str                     | output base name
# out
#   bcf      file  =  ${outbase}.bcf | output alignment
# run
#   time   = 10
#   image  = djhshih/seqkit:0.1

bcftools merge -m both -o ${outbase}.bcf ${bcf}
