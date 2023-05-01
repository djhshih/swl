#? Call mutations on read alignment
#
# in  bam      file                    | input bam
# in  ref      file                    | reference sequence
# in  ref_fai  file                    | reference index
# in  outbase  str                     | output base name
# out bcf      file  =  ${outbase}.bcf | output alignment
# run image  = djhshih/seqkit:0.1

bcftools mpileup -Ou -f ${ref} ${bam} | 
	bcftools call -mv -Ob -o ${bcf}
