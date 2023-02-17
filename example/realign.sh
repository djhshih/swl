#' Align paired-end sequencing reads
#' @in  fastq1   file    read 1
#' @in  fastq2   file    read 2
#' @in  ref      file    reference sequence
#' @in  bam      str     output file name
#' @out bam      file
#' @out bai      file
#' @run cpu      int     [default: 2]
#' @run memory   int     [default: 8G]

bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > aligned.bam

samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${bam}
samtools index ${bam} ${bai}
