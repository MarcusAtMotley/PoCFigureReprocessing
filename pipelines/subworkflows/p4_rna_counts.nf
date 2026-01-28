/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P4: RNA Gene Counting Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Counts gene expression from aligned RNA BAMs with FeatureCounts.

    Input:  Aligned BAM + BAI (from ALIGN_RNA)
    Output: Gene count matrix

    Flow: BAM → FeatureCounts → Counts CSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SUBREAD_FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'

workflow P4_RNA_COUNTS {

    take:
    ch_bam  // channel: [ val(meta), bam ]
    ch_bai  // channel: [ val(meta), bai ]
    ch_gtf  // channel: [ val(meta), gtf ]

    main:
    ch_versions = Channel.empty()

    //
    // Prepare input for FeatureCounts
    // Expects: [ meta, bam, gtf ]
    // Note: combine with full ch_gtf (value channel) to enable streaming
    //
    ch_featurecounts_input = ch_bam
        .combine(ch_gtf)
        .map { meta, bam, meta2, gtf -> [ meta, bam, gtf ] }

    //
    // MODULE: FeatureCounts for gene quantification
    //
    SUBREAD_FEATURECOUNTS (
        ch_featurecounts_input
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    emit:
    counts   = SUBREAD_FEATURECOUNTS.out.counts   // channel: [ val(meta), counts.tsv ]
    summary  = SUBREAD_FEATURECOUNTS.out.summary  // channel: [ val(meta), summary ]
    versions = ch_versions                         // channel: [ versions.yml ]
}
