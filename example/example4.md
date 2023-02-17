Tasks are annotated bash scripts.

`align.sh`
```{bash}
#' doc Align paired-end sequencing reads
#' in  fastq1   file                    | read 1
#' in  fastq2   file                    | read 2
#' in  ref      file                    | reference sequence
#' in  outbase  str                     | output base name
#' out bam      file  =  ${outbase}.bam | output alignment
#' run cpu      int   = 2
#' run image    size  = djhshih/seqkit:0.1

bwa mem -t ${cpu} ${ref} ${fastq1} ${fastq2} | samtools view -b - > ${outbase}.bam
```

`sort.sh`
```{bash}
#' doc Sort alignment by coordinates and index
#' in  bam      file                    | input bam
#' in  ref      file                    | reference sequence
#' in  outbase  str                     | output base name
#' out bam      file  =  ${outbase}.bam | output alignment
#' out bai      file  =  ${outbase}.bai | output alignment index
#' run cpu    = 2
#' run memory = 8G
#' run image  = djhshih/seqkit:0.1

samtools sort -@ ${cpu} -m ${memory / cpu} aligned.bam ${outbase}.bam
samtools index ${outbase}.bam ${outbase}.bai
```

Tasks output records.
From above, we see that task `align` takes an input record of 
`fastq1` and `fastq2`, and it outputs a record of `bam`, and 
task `sort` takes a record of `bam` and outputs a record of `bam` and `bai`. 

We can import the tasks into a workflow by

`align.workflow`
```
align = import "align.sh"
sort  = import "sort.sh"
```

We can then write a workflow as
```
align |> sort
```
which is syntactic sugar for
```
a = align argv
s = sort { bam = a.bam }

{ bam = s.bam, bai = s.bai }
```
where `argv` is a record of the workflow input. Comparing this two blocks of code, we also note that

- Each task is separated from its input record by a space.
- Records are constructed using `{ }`.
- The first task is automatically given the workflow input `argv`.
- Each subsequent task produce a record. Attributes in later records take precedence.
- The output of a workflow is just the resulting record on the last line.

Attributes in the input record may be qualified to resolve ambiguity:
```
{ align = { fastq1 = "1.fq", fastq2 = "2.fq" } }
```

Additionally, we can merge records by
```
a & s
```
The record on the right is always right (i.e. it takes precedence when both records have the same attribute).

This workflow can also be imported into another workflow similarly as above:
```
align = import "align.workflow"
```

Documentation for the workflow can be generated using the documentation of the tasks.
