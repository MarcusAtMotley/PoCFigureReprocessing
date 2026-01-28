/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P3: CNV Calling Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls copy number variants from aligned DNA BAMs using CNVpytor,
    with optional intergenic filtering.

    Input:  Aligned BAM + BAI (from ALIGN_DNA)
    Output: CNV calls in CSV format

    Flow: BAM → CNVpytor → Intergenic Filter → CNV CSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CNVPYTOR_CALL     } from '../modules/local/cnvpytor/main'
include { INTERGENIC_FILTER } from '../modules/local/intergenic_filter/main'

workflow P3_CNV {

    take:
    ch_bam            // channel: [ val(meta), bam ]
    ch_bai            // channel: [ val(meta), bai ]
    ch_fasta          // channel: [ val(meta), fasta ]
    ch_fai            // channel: [ val(meta), fai ]
    ch_intergenic_bed // channel: path(intergenic.bed) or empty

    main:
    ch_versions = Channel.empty()

    //
    // Prepare BAM channel with index for CNVpytor
    //
    ch_bam_for_cnv = ch_bam.join(ch_bai)

    //
    // MODULE: CNVpytor CNV calling
    //
    CNVPYTOR_CALL (
        ch_bam_for_cnv,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(CNVPYTOR_CALL.out.versions.first())

    //
    // MODULE: Intergenic filter (optional)
    // Only run if intergenic BED file is provided
    //
    if (ch_intergenic_bed) {
        INTERGENIC_FILTER (
            CNVPYTOR_CALL.out.calls,
            ch_intergenic_bed
        )
        ch_versions = ch_versions.mix(INTERGENIC_FILTER.out.versions.first())
        ch_cnv_output = INTERGENIC_FILTER.out.filtered_calls
    } else {
        // No filtering, use CNVpytor output directly
        ch_cnv_output = CNVPYTOR_CALL.out.calls
    }

    emit:
    pytor    = CNVPYTOR_CALL.out.pytor  // channel: [ val(meta), pytor ]
    cnv      = ch_cnv_output            // channel: [ val(meta), cnv_csv ]
    versions = ch_versions               // channel: [ versions.yml ]
}
