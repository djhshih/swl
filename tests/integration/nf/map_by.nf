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
    xs_ch = Channel.fromList(params.xs)

    call_variant_ch = xs_ch
        .map{ x -> tuple(
            x.outbase,
            file(x.fastq1),
            file(x.fastq2),
            file(x.ref),
            file(x.ref_amb),
            file(x.ref_ann),
            file(x.ref_bwt),
            file(x.ref_fai),
            file(x.ref_pac),
            file(x.ref_sa)
        ) }
        .groupTuple()

    CALL_VARIANT(call_variant_ch)


    bai = CALL_VARIANT.out.bai
    bam = CALL_VARIANT.out.bam
    bcf = CALL_VARIANT.out.bcf
}
