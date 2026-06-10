process ALIGN {

    input:
    path fastq1
    path fastq2
    val outbase
    path ref
    path ref_amb
    path ref_ann
    path ref_bwt
    path ref_pac
    path ref_sa

    output:
    path "${outbase}.bam", emit: bam

    script:
    """
    cat "${fastq1}" "${fastq2}" "${ref}" "${ref_amb}" "${ref_ann}" "${ref_bwt}" "${ref_pac}" "${ref_sa}" > "${outbase}.bam"
    """

}

process SORT {

    input:
    path bam
    val outbase

    output:
    path "${outbase}.bai", emit: bai
    path "${outbase}.bam", emit: bam

    script:
    """
    cp "${bam}" "${outbase}.bam"
    cp "${bam}" "${outbase}.bai"
    """

}

process CALL {

    input:
    path bam
    val outbase
    path ref
    path ref_fai

    output:
    path "${outbase}.bcf", emit: bcf

    script:
    """
    cat "${bam}" "${ref}" "${ref_fai}" > "${outbase}.bcf"
    """

}

workflow {
    fastq1_ch = Channel.fromPath(params.fastq1, checkIfExists: true)
    fastq2_ch = Channel.fromPath(params.fastq2, checkIfExists: true)
    outbase_ch = Channel.value(params.outbase)
    ref_ch = Channel.fromPath(params.ref, checkIfExists: true)
    ref_amb_ch = Channel.fromPath(params.ref_amb, checkIfExists: true)
    ref_ann_ch = Channel.fromPath(params.ref_ann, checkIfExists: true)
    ref_bwt_ch = Channel.fromPath(params.ref_bwt, checkIfExists: true)
    ref_fai_ch = Channel.fromPath(params.ref_fai, checkIfExists: true)
    ref_pac_ch = Channel.fromPath(params.ref_pac, checkIfExists: true)
    ref_sa_ch = Channel.fromPath(params.ref_sa, checkIfExists: true)

    ALIGN(fastq1_ch, fastq2_ch, outbase_ch, ref_ch, ref_amb_ch, ref_ann_ch, ref_bwt_ch, ref_pac_ch, ref_sa_ch)
    SORT(ALIGN.out.bam, outbase_ch)
    CALL(SORT.out.bam, outbase_ch, ref_ch, ref_fai_ch)

    bai = SORT.out.bai
    bam = SORT.out.bam
    bcf = CALL.out.bcf
}
