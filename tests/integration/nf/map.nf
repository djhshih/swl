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

    fastq1 = xs_ch.map{ x -> file(x.fastq1) }
    fastq2 = xs_ch.map{ x -> file(x.fastq2) }
    outbase = xs_ch.map{ x -> x.outbase }
    ref = xs_ch.map{ x -> file(x.ref) }
    ref_amb = xs_ch.map{ x -> file(x.ref_amb) }
    ref_ann = xs_ch.map{ x -> file(x.ref_ann) }
    ref_bwt = xs_ch.map{ x -> file(x.ref_bwt) }
    ref_fai = xs_ch.map{ x -> file(x.ref_fai) }
    ref_pac = xs_ch.map{ x -> file(x.ref_pac) }
    ref_sa = xs_ch.map{ x -> file(x.ref_sa) }

    ALIGN(fastq1, fastq2, outbase, ref, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa)
    SORT(ALIGN.out.bam, outbase)
    CALL(SORT.out.bam, outbase, ref, ref_fai)


    bai = SORT.out.bai.toList()
    bam = SORT.out.bam.toList()
    bcf = CALL.out.bcf.toList()
}
