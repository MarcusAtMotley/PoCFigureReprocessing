/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DNA Alignment Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared Biscuit alignment for all DNA pipelines (P1, P2, P3).
    Align once, use for SNP calling, methylation, and CNV analysis.

    Input:  PE or SE FASTQs (WGS, WGEM, mTNA-DNA, HairyTNA-DNA)
    Output: Sorted BAM + BAI for downstream analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BISCUIT_ALIGN } from '../modules/nf-core/biscuit/align/main'

workflow ALIGN_DNA {

    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_fasta         // channel: [ val(meta), fasta ]
    ch_biscuit_index // channel: [ val(meta), index_dir ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Biscuit alignment
    // Bisulfite-aware alignment for DNA reads
    //
    BISCUIT_ALIGN (
        ch_reads,
        ch_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_ALIGN.out.versions.first())

    emit:
    bam      = BISCUIT_ALIGN.out.bam     // channel: [ val(meta), bam ]
    bai      = BISCUIT_ALIGN.out.bai     // channel: [ val(meta), bai ]
    versions = ch_versions                // channel: [ versions.yml ]
}
