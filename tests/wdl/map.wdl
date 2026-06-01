version 1.1

version 1.1

task align {

    input {
        File fastq1
        File fastq2
        String outbase
        File ref
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
    }

    command <<<
        bwa mem -t ~{cpu} ~{ref} ~{fastq1} ~{fastq2} | samtools view -b - > ~{outbase}.bam
        
    >>>

    output {
        File bam = "~{outbase}.bam"
    }

    requirements {
        cpu: 2
        container: "djhshih/seqkit:0.1"
        time_minutes: 3150
    }

}

task sort {

    input {
        File bam
        String outbase
    }

    command <<<
        samtools sort -@ ~{cpu} -m ~{memory / cpu} aligned.bam ~{outbase}.bam
        samtools index ~{outbase}.bam ~{outbase}.bai
        
    >>>

    output {
        File bai = "~{outbase}.bai"
        File bam = "~{outbase}.bam"
    }

    requirements {
        cpu: 2
        container: "djhshih/seqkit:0.1"
        memory: "8192 MB"
    }

}

task call {

    input {
        File bam
        String outbase
        File ref
        File ref_fai
    }

    command <<<
        bcftools mpileup -Ou -f ~{ref} ~{bam} | 
        	bcftools call -mv -Ob -o ~{outbase}.bcf
        bcftools index ~{outbase}.bcf
        
    >>>

    output {
        File bcf = "~{outbase}.bcf"
    }

    requirements {
        container: "djhshih/seqkit:0.1"
        time_minutes: 30
    }

}

workflow main_call_variant {

    input {
        File fastq1
        File fastq2
        String outbase
        File ref
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_fai
        File ref_pac
        File ref_sa
    }

    call align {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            outbase = outbase,
            ref = ref,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
    }

    call sort {
        input:
            outbase = outbase,
            bam = align.bam,
    }

    call call {
        input:
            outbase = outbase,
            ref = ref,
            ref_fai = ref_fai,
            bam = sort.bam,
    }

    output {
        File bai = sort.bai
        File bam = sort.bam
        File bcf = call.bcf
    }
}

workflow main {

    input {
        String xs
    }

    scatter (call_variant_i in range(length(xs))) {
        call call_variant {
            input:
                fastq1 = fastq1,
                fastq2 = fastq2,
                outbase = outbase,
                ref = ref,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_bwt = ref_bwt,
                ref_fai = ref_fai,
                ref_pac = ref_pac,
                ref_sa = ref_sa,
        }
    }

    output {
        Array[File] bai = call_variant.bai
        Array[File] bam = call_variant.bam
        Array[File] bcf = call_variant.bcf
    }
}
