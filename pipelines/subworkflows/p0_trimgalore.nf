/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P0: TrimGalore Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Trims single-analyte raw FASTQs using TrimGalore.
    Only needed for samples from s3://motleybio-medgenome/ (raw data).

    Input:  Raw paired-end FASTQs
    Output: Trimmed FASTQs ready for downstream pipelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'

workflow P0_TRIMGALORE {

    take:
    ch_reads  // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run TrimGalore for quality and adapter trimming
    //
    TRIMGALORE ( ch_reads )
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    emit:
    reads       = TRIMGALORE.out.reads    // channel: [ val(meta), [ trimmed_reads ] ]
    trim_log    = TRIMGALORE.out.log      // channel: [ val(meta), log ]
    html        = TRIMGALORE.out.html     // channel: [ val(meta), html ]
    zip         = TRIMGALORE.out.zip      // channel: [ val(meta), zip ]
    versions    = ch_versions             // channel: [ versions.yml ]
}
