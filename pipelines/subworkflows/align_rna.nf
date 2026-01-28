/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RNA Alignment Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared STAR alignment for all RNA pipelines (P4, P5).
    Align once, use for both FeatureCounts and SNP calling.

    Input:  PE or SE FASTQs (RNA-Seq, mTNA-RNA, HairyTNA-RNA)
    Output: Sorted BAM + BAI for downstream analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAR_ALIGN     } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

workflow ALIGN_RNA {

    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_star_index // channel: [ val(meta), index_dir ]
    ch_gtf        // channel: [ val(meta), gtf ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: STAR alignment
    // Splice-aware alignment for RNA reads
    //
    STAR_ALIGN (
        ch_reads,
        ch_star_index,
        ch_gtf,
        false,  // star_ignore_sjdbgtf
        '',     // seq_platform
        ''      // seq_center
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // MODULE: Index the sorted BAM
    //
    SAMTOOLS_INDEX (
        STAR_ALIGN.out.bam_sorted
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = STAR_ALIGN.out.bam_sorted    // channel: [ val(meta), bam ]
    bai      = SAMTOOLS_INDEX.out.bai       // channel: [ val(meta), bai ]
    log      = STAR_ALIGN.out.log_final     // channel: [ val(meta), log ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
