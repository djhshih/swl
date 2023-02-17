#' doc Align paired-end sequencing reads
#' in  fastq1   file                    | read 1
#' in  fastq2   file                    | read 2
#' in  ref      file                    | reference sequence
#' in  outbase  str                     | output base name
#' out bam      file  =  {outbase}.bam  | output alignment
#' out bai      file  =  {outbase}.bai  | output alignment index
#' run cpu    = 2
#' run memory = 8G
#' run image  = djhshih/seqkit:0.1

bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > aligned.bam

samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${bam}
samtools index ${bam} ${bai}
