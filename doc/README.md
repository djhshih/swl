Tasks are annotated bash scripts.

`align.sh`
```{bash}
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
# run image    str   = djhshih/seqkit:0.1

bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam
```

`sort.sh`
```{bash}
#? Sort alignment by coordinates and index
#
# in  bam      file                    | input bam
# in  outbase  str                     | output base name
# out bam      file  =  ${outbase}.bam | output alignment
# out bai      file  =  ${outbase}.bai | output alignment index
# run cpu    = 2
# run memory = 8G
# run image  = djhshih/seqkit:0.1

samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam
samtools index ${outbase}.bam ${outbase}.bai
```

`call.sh`
```{bash}
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
```

Tasks output records.
From above, we see that task `align` takes an input record of 
`fastq1` and `fastq2`, and it outputs a record of `bam`, and 
task `sort` takes a record of `bam` and outputs a record of `bam` and `bai`. 

We can import the tasks into a workflow by

`mutation.wf`
```
align = import "align.sh"
sort  = import "sort.sh"
call  = import "call.sh"
```

We can then write the workflow in `mutation.wf` simply with
```
align |> sort |> call
```
which is syntactic sugar for
```
\argv ->
    a = align argv
    s = sort ( argv & a )
    c = call ( argv & a & s )
    a & s & c
```
where `argv` is a record of the workflow input, and `\argv -> block` is a function.
The result of a block of code is simply given by the final line.

Comparing this two blocks of code, we also note that
- Each task is separated from its input record by a space.
- Records are constructed using `{ }`.
- The first task is run with the workflow input `argv`, producing a record that is then given to subsequent tasks.
- Each task produce a record using the workflow input merged with output records from the preceding tasks.
- The output of a workflow is result on the final line.

As seen above, we merge records by
```
s & c
```
The record on the right is always right (i.e. it takes precedence when both records have the same attribute).

Attributes in the input `argv` record may be qualified to resolve ambiguity:
```
{ align: { fastq1: "1.fq", fastq2: "2.fq" } }
```

Alternatively, we can also write the workflow with explicit passing of arguments by
```
\argv ->
    a = align argv
    s = sort { bam: a.bam, outbase: argv.outbase }
    c = call { bam: s.bam, ref: argv.ref, ref_fai: argv.ref_fai, outbase: argv.outbase }

    { bam: s.bam, bai: s.bai, bcf: c.bcf }
```

A workflow can be imported into another workflow similarly as above:
```
mutation = import "mutation.wf"
```

Documentation for the workflow can be generated using the documentation of the tasks.

### Partial function application

When a function is applied to a record that does not have the full set of attributes,
another function is returned.
This function can then be applied to another record.
```
align = import "align.sh"

align_hg38 = align {
  ref: "hg38.fa",
  ref_ann: "hg38.fa.ann",
  ref_bwt: "hg38.fa.bwt",
  ref_pac: "hg38.fa.pac",
  ref_sa: "hg38.fa.sa",
}

align_hg38 { fastq1: "sample1.r1.fq", fastq2: "sample.r2.fq" }
```
