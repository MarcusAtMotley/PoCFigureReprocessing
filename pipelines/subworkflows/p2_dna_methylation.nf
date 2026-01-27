/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P2: DNA Methylation Calling Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls methylation from DNA channel reads using Biscuit alignment and pileup.

    Input:  PE or SE FASTQs (WGEM, mTNA-DNA, HairyTNA-DNA)
    Output: Methylation VCF

    Flow: FASTQ → Biscuit Align → Biscuit Pileup → Methylation VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BISCUIT_ALIGN  } from '../modules/nf-core/biscuit/align/main'
include { BISCUIT_PILEUP } from '../modules/nf-core/biscuit/pileup/main'

workflow P2_DNA_METH {

    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_fasta         // channel: [ val(meta), fasta ]
    ch_biscuit_index // channel: [ val(meta), index_dir ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Biscuit alignment
    // Aligns reads to reference genome (bisulfite-aware)
    //
    BISCUIT_ALIGN (
        ch_reads,
        ch_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_ALIGN.out.versions.first())

    //
    // Prepare BAM channel for Biscuit pileup
    // Pileup expects: [ meta, normal_bams, normal_bais, tumor_bam, tumor_bai ]
    // For single-sample methylation calling, pass empty tumor channels
    //
    ch_bam_for_pileup = BISCUIT_ALIGN.out.bam
        .join(BISCUIT_ALIGN.out.bai)
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
    bam      = BISCUIT_ALIGN.out.bam     // channel: [ val(meta), bam ]
    bai      = BISCUIT_ALIGN.out.bai     // channel: [ val(meta), bai ]
    vcf      = BISCUIT_PILEUP.out.vcf    // channel: [ val(meta), vcf.gz ]
    versions = ch_versions                // channel: [ versions.yml ]
}
