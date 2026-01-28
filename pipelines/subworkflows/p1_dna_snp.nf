/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P1: DNA SNP Calling Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls SNPs from aligned DNA BAMs using Revelio and LoFreq.

    Input:  Aligned BAM + BAI (from ALIGN_DNA)
    Output: VCF with SNP calls

    Flow: BAM → Revelio → LoFreq → VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REVELIO              } from '../modules/local/revelio/main'
include { LOFREQ_CALLPARALLEL  } from '../modules/nf-core/lofreq/callparallel/main'

workflow P1_DNA_SNP {

    take:
    ch_bam   // channel: [ val(meta), bam ]
    ch_bai   // channel: [ val(meta), bai ]
    ch_fasta // channel: [ val(meta), fasta ]
    ch_fai   // channel: [ val(meta), fai ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Revelio - mask bisulfite conversions
    // Transforms BAM to allow standard variant calling on bisulfite data
    //
    ch_bam_for_revelio = ch_bam.join(ch_bai)

    REVELIO (
        ch_bam_for_revelio,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(REVELIO.out.versions.first())

    //
    // Prepare BAM channel with index for LoFreq
    // LoFreq expects: [ meta, bam, bai, intervals ]
    //
    ch_bam_for_lofreq = REVELIO.out.bam
        .join(REVELIO.out.bai)
        .map { meta, bam, bai -> [ meta, bam, bai, [] ] }

    //
    // MODULE: LoFreq variant calling
    //
    LOFREQ_CALLPARALLEL (
        ch_bam_for_lofreq,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(LOFREQ_CALLPARALLEL.out.versions.first())

    emit:
    bam      = REVELIO.out.bam                   // channel: [ val(meta), bam ] (revelio-masked)
    bai      = REVELIO.out.bai                   // channel: [ val(meta), bai ]
    vcf      = LOFREQ_CALLPARALLEL.out.vcf       // channel: [ val(meta), vcf.gz ]
    tbi      = LOFREQ_CALLPARALLEL.out.tbi       // channel: [ val(meta), vcf.gz.tbi ]
    versions = ch_versions                        // channel: [ versions.yml ]
}
