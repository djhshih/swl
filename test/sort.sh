#? Sort alignment by coordinates and index
# in  bam      file                    | input bam
# in  outbase  str                     | output base name
# out bam      file  =  ${outbase}.bam | output alignment
# out bai      file  =  ${outbase}.bai | output alignment index
# run cpu    = 2
# run memory = 8G
# run image  = djhshih/seqkit:0.1

samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam
samtools index ${outbase}.bam ${outbase}.bai
