/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P4: RNA Gene Counting Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Aligns RNA reads with STAR and counts gene expression with FeatureCounts.

    Input:  PE or SE FASTQs (RNA-Seq, mTNA-RNA, HairyTNA-RNA)
    Output: Gene count matrix

    Flow: FASTQ → STAR Align → FeatureCounts → Counts CSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAR_ALIGN            } from '../modules/nf-core/star/align/main'
include { SUBREAD_FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'

workflow P4_RNA_COUNTS {

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
    // Prepare BAM channel for FeatureCounts
    // Use sorted BAM output from STAR
    //
    ch_bam_for_counts = STAR_ALIGN.out.bam_sorted
        .map { meta, bam -> [ meta, bam, ch_gtf.map{ it[1] }.first() ] }

    // Join BAM with GTF for featurecounts input
    ch_featurecounts_input = STAR_ALIGN.out.bam_sorted
        .combine(ch_gtf.map { meta, gtf -> gtf })
        .map { meta, bam, gtf -> [ meta, bam, gtf ] }

    //
    // MODULE: FeatureCounts for gene quantification
    //
    SUBREAD_FEATURECOUNTS (
        ch_featurecounts_input
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    emit:
    bam         = STAR_ALIGN.out.bam_sorted           // channel: [ val(meta), bam ]
    counts      = SUBREAD_FEATURECOUNTS.out.counts    // channel: [ val(meta), counts.tsv ]
    summary     = SUBREAD_FEATURECOUNTS.out.summary   // channel: [ val(meta), summary ]
    star_log    = STAR_ALIGN.out.log_final            // channel: [ val(meta), log ]
    versions    = ch_versions                          // channel: [ versions.yml ]
}
