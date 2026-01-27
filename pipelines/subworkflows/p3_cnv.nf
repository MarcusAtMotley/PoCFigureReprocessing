/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P3: CNV Calling Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls copy number variants from DNA channel reads using Biscuit alignment
    and CNVpytor, with optional intergenic filtering.

    Input:  PE or SE FASTQs (WGS, mTNA-DNA, HairyTNA-DNA)
    Output: CNV calls in CSV format

    Flow: FASTQ → Biscuit Align → CNVpytor → Intergenic Filter → CNV CSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BISCUIT_ALIGN     } from '../modules/nf-core/biscuit/align/main'
include { CNVPYTOR_CALL     } from '../modules/local/cnvpytor/main'
include { INTERGENIC_FILTER } from '../modules/local/intergenic_filter/main'

workflow P3_CNV {

    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_fasta         // channel: [ val(meta), fasta ]
    ch_fai           // channel: [ val(meta), fai ]
    ch_biscuit_index // channel: [ val(meta), index_dir ]
    ch_intergenic_bed // channel: path(intergenic.bed) or empty

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Biscuit alignment
    // Aligns reads to reference genome
    //
    BISCUIT_ALIGN (
        ch_reads,
        ch_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_ALIGN.out.versions.first())

    //
    // Prepare BAM channel with index for CNVpytor
    //
    ch_bam_for_cnv = BISCUIT_ALIGN.out.bam
        .join(BISCUIT_ALIGN.out.bai)

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
        // No filtering, convert TSV to CSV directly
        ch_cnv_output = CNVPYTOR_CALL.out.calls.map { meta, tsv ->
            [ meta, tsv ]
        }
    }

    emit:
    bam      = BISCUIT_ALIGN.out.bam    // channel: [ val(meta), bam ]
    bai      = BISCUIT_ALIGN.out.bai    // channel: [ val(meta), bai ]
    pytor    = CNVPYTOR_CALL.out.pytor  // channel: [ val(meta), pytor ]
    cnv      = ch_cnv_output            // channel: [ val(meta), cnv_csv ]
    versions = ch_versions               // channel: [ versions.yml ]
}
