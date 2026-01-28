/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P2: DNA Methylation Calling Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls methylation from aligned DNA BAMs using Biscuit Pileup.

    Input:  Aligned BAM + BAI (from ALIGN_DNA)
    Output: Methylation VCF

    Flow: BAM → Biscuit Pileup → Methylation VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BISCUIT_PILEUP } from '../modules/nf-core/biscuit/pileup/main'

workflow P2_DNA_METH {

    take:
    ch_bam           // channel: [ val(meta), bam ]
    ch_bai           // channel: [ val(meta), bai ]
    ch_fasta         // channel: [ val(meta), fasta ]
    ch_biscuit_index // channel: [ val(meta), index_dir ]

    main:
    ch_versions = Channel.empty()

    //
    // Prepare BAM channel for Biscuit pileup
    // Pileup expects: [ meta, normal_bams, normal_bais, tumor_bam, tumor_bai ]
    // For single-sample methylation calling, pass empty tumor channels
    //
    ch_bam_for_pileup = ch_bam
        .join(ch_bai)
        .map { meta, bam, bai -> [ meta, bam, bai, [], [] ] }

    //
    // MODULE: Biscuit pileup for methylation calling
    //
    BISCUIT_PILEUP (
        ch_bam_for_pileup,
        ch_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_PILEUP.out.versions.first())

    emit:
    vcf      = BISCUIT_PILEUP.out.vcf    // channel: [ val(meta), vcf.gz ]
    versions = ch_versions                // channel: [ versions.yml ]
}
