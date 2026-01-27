/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    P5: RNA SNP Calling Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls SNPs from RNA reads using STAR alignment and LoFreq.

    Input:  PE or SE FASTQs (RNA-Seq, mTNA-RNA, HairyTNA-RNA)
    Output: VCF with RNA SNP calls

    Flow: FASTQ → STAR Align → Index → LoFreq → VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAR_ALIGN           } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX       } from '../modules/nf-core/samtools/index/main'
include { LOFREQ_CALLPARALLEL  } from '../modules/nf-core/lofreq/callparallel/main'

workflow P5_RNA_SNP {

    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_star_index // channel: [ val(meta), index_dir ]
    ch_gtf        // channel: [ val(meta), gtf ]
    ch_fasta      // channel: [ val(meta), fasta ]
    ch_fai        // channel: [ val(meta), fai ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: STAR alignment
    // Splice-aware alignment for RNA reads
    // Using two-pass mode for better splice junction detection
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

    //
    // Prepare BAM channel with index for LoFreq
    // LoFreq expects: [ meta, bam, bai, intervals ]
    //
    ch_bam_for_lofreq = STAR_ALIGN.out.bam_sorted
        .join(SAMTOOLS_INDEX.out.bai)
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
    bam         = STAR_ALIGN.out.bam_sorted       // channel: [ val(meta), bam ]
    bai         = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), bai ]
    vcf         = LOFREQ_CALLPARALLEL.out.vcf     // channel: [ val(meta), vcf.gz ]
    tbi         = LOFREQ_CALLPARALLEL.out.tbi     // channel: [ val(meta), vcf.gz.tbi ]
    star_log    = STAR_ALIGN.out.log_final        // channel: [ val(meta), log ]
    versions    = ch_versions                      // channel: [ versions.yml ]
}
