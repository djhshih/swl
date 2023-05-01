#? Align paired-end sequencing reads
#
# in  fastq1   file                    | read 1
# in  fastq2   file                    | read 2
# in  ref      file                    | reference sequence
# in  ref_amb  file                    | reference bwa index file
# in  ref_ann  file                    | reference bwa index file
# in  ref_bwt  file                    | reference bwa index file
# in  ref_pac  file                    | reference bwa index file
# in  ref_sa   file                    | reference bwa index file
# in  outbase  str                     | output base name
# out bam      file  =  ${outbase}.bam | output alignment
# run cpu      int   = 2
# run image    size  = djhshih/seqkit:0.1

bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam
