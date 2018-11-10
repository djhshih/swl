bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > aligned.bam

samtools sort -@ ${cpu} -m ${memory/cpu} aligned.bam ${bam}
samtools index ${bam} ${bai}
