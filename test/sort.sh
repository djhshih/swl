#?  Sort alignment by coordinates and index
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
