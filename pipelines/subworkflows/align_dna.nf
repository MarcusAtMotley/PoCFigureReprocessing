/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DNA Alignment Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared Biscuit alignment for all DNA pipelines (P1, P2, P3).
    Align once, use for SNP calling, methylation, and CNV analysis.

    Input:  PE or SE FASTQs (WGS, WGEM, mTNA-DNA, HairyTNA-DNA)
    Output: Sorted BAM + BAI for downstream analysis

    Uses scatter-gather for large files:
    1. Split FASTQs into chunks (configurable, default 10M read pairs)
    2. Align each chunk in parallel
    3. Merge BAMs back together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BISCUIT_ALIGN  } from '../modules/nf-core/biscuit/align/main'
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

workflow ALIGN_DNA {

    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_fasta         // channel: [ val(meta), fasta ]
    ch_biscuit_index // channel: [ val(meta), index_dir ]

    main:
    ch_versions = Channel.empty()

    // Get chunk size from params (default 10M read pairs)
    // Set to 0 to disable scatter-gather
    def chunk_size = params.biscuit_chunk_size ?: 10000000

    if (chunk_size > 0) {
        //
        // SCATTER-GATHER MODE: Split large FASTQs for parallel alignment
        //

        // Branch by single_end status since splitFastq syntax differs
        ch_reads
            .branch { meta, reads ->
                pe: !meta.single_end && reads.size() > 1
                se: meta.single_end || reads.size() == 1
            }
            .set { ch_branched }

        // Split paired-end reads
        // elem: [1] tells splitFastq which tuple element contains the files
        ch_pe_split = ch_branched.pe
            .map { meta, reads -> [ meta, reads[0], reads[1] ] }
            .splitFastq(by: chunk_size, pe: true, elem: [1, 2], file: true)
            .map { meta, r1, r2 -> [ meta, [r1, r2] ] }

        // Split single-end reads
        ch_se_split = ch_branched.se
            .map { meta, reads ->
                def r1 = reads instanceof List ? reads[0] : reads
                [ meta, r1 ]
            }
            .splitFastq(by: chunk_size, elem: 1, file: true)
            .map { meta, r1 -> [ meta, [r1] ] }

        // Combine split channels
        ch_reads_chunked = ch_pe_split.mix(ch_se_split)

        //
        // MODULE: Biscuit alignment on each chunk
        //
        BISCUIT_ALIGN (
            ch_reads_chunked,
            ch_fasta,
            ch_biscuit_index
        )
        ch_versions = ch_versions.mix(BISCUIT_ALIGN.out.versions.first())

        //
        // GATHER: Group BAMs by sample ID and merge
        //
        ch_bams_grouped = BISCUIT_ALIGN.out.bam
            .map { meta, bam ->
                // Create a clean key for grouping (strip split index from meta)
                def key = meta.id
                [ key, meta, bam ]
            }
            .groupTuple(by: 0)
            .map { key, metas, bams ->
                // Use first meta, flatten bams list
                [ metas[0], bams.flatten() ]
            }

        // Merge all chunks (works for single chunk too - just passes through)
        ch_fasta_for_merge = ch_fasta.map { meta, fa -> [ [:], fa ] }
        ch_empty = Channel.value([ [:], [] ])

        SAMTOOLS_MERGE (
            ch_bams_grouped,
            ch_fasta_for_merge.first(),
            ch_empty,
            ch_empty
        )

        // Index merged BAMs
        SAMTOOLS_INDEX ( SAMTOOLS_MERGE.out.bam )

        ch_final_bam = SAMTOOLS_MERGE.out.bam
        ch_final_bai = SAMTOOLS_INDEX.out.bai

    } else {
        //
        // SIMPLE MODE: No scatter-gather, direct alignment
        //
        BISCUIT_ALIGN (
            ch_reads,
            ch_fasta,
            ch_biscuit_index
        )
        ch_versions = ch_versions.mix(BISCUIT_ALIGN.out.versions.first())

        ch_final_bam = BISCUIT_ALIGN.out.bam
        ch_final_bai = BISCUIT_ALIGN.out.bai
    }

    emit:
    bam      = ch_final_bam   // channel: [ val(meta), bam ]
    bai      = ch_final_bai   // channel: [ val(meta), bai ]
    versions = ch_versions    // channel: [ versions.yml ]
}
