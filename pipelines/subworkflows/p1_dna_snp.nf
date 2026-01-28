/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P1: DNA SNP Calling Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls SNPs from DNA channel reads using Biscuit alignment, Revelio bisulfite
    masking, and LoFreq variant calling.

    Input:  PE or SE FASTQs (WGS, mTNA-DNA, HairyTNA-DNA)
    Output: VCF with SNP calls

    Flow: FASTQ → Biscuit Align → Revelio → LoFreq → VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BISCUIT_ALIGN        } from '../modules/nf-core/biscuit/align/main'
include { REVELIO              } from '../modules/local/revelio/main'
include { LOFREQ_CALLPARALLEL  } from '../modules/nf-core/lofreq/callparallel/main'

workflow P1_DNA_SNP {

    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_fasta         // channel: [ val(meta), fasta ]
    ch_fai           // channel: [ val(meta), fai ]
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
    // MODULE: Revelio - mask bisulfite conversions
    // Transforms BAM to allow standard variant calling on bisulfite data
    //
    ch_bam_for_revelio = BISCUIT_ALIGN.out.bam
        .join(BISCUIT_ALIGN.out.bai)

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
        .map { meta, bam, bai -> [ meta, bam, bai, [] ] }  // Empty intervals = call on all regions

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
