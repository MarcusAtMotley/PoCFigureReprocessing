/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DNA Alignment Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared Biscuit alignment for all DNA pipelines (P1, P2, P3).
    Align once, use for SNP calling, methylation, and CNV analysis.

    Input:  PE or SE FASTQs (WGS, WGEM, mTNA-DNA, HairyTNA-DNA)
    Output: Sorted BAM + BAI for downstream analysis

    Uses scatter-gather for large files:
    1. Split FASTQs into chunks (default 10M read pairs)
    2. Align each chunk in parallel
    3. Cat + sort + index merged BAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BISCUIT_ALIGN  } from '../modules/nf-core/biscuit/align/main'
include { SAMTOOLS_SORT  } from '../modules/nf-core/samtools/sort/main'

workflow ALIGN_DNA {

    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_fasta         // channel: [ val(meta), fasta ]
    ch_biscuit_index // channel: [ val(meta), index_dir ]

    main:
    ch_versions = Channel.empty()

    // Get chunk size from params (default 10M read pairs)
    def chunk_size = params.biscuit_chunk_size ?: 10000000

    //
    // SCATTER: Split FASTQs into chunks for parallel alignment
    // Branch by single_end status since splitFastq syntax differs
    //
    ch_reads
        .branch { meta, reads ->
            pe: !meta.single_end
            se: meta.single_end
        }
        .set { ch_branched }

    // Split paired-end reads (keeps R1/R2 synchronized)
    ch_pe_split = ch_branched.pe
        .splitFastq(by: chunk_size, pe: true, file: true, compress: true)

    // Split single-end reads
    ch_se_split = ch_branched.se
        .splitFastq(by: chunk_size, file: true, compress: true)

    // Combine back into single channel
    ch_reads_chunked = ch_pe_split.mix(ch_se_split)

    //
    // MODULE: Biscuit alignment on each chunk
    // Note: BISCUIT_ALIGN already sorts each chunk individually
    //
    BISCUIT_ALIGN (
        ch_reads_chunked,
        ch_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_ALIGN.out.versions.first())

    //
    // GATHER: Group BAMs by sample ID, then cat + sort + index
    // SAMTOOLS_SORT does: samtools cat | samtools sort --write-index
    //
    ch_bams_to_merge = BISCUIT_ALIGN.out.bam
        .map { meta, bam ->
            // Strip any chunk suffix from meta for grouping
            def clean_meta = meta.clone()
            clean_meta.remove('split')  // splitFastq adds 'split' to meta
            [ clean_meta.id, clean_meta, bam ]
        }
        .groupTuple(by: 0)
        .map { id, metas, bams ->
            // Use first meta (they're all the same except for split info)
            [ metas[0], bams.flatten() ]
        }

    // Reference for sorting (optional but good for CRAM compatibility)
    ch_fasta_for_sort = ch_fasta.map { meta, fasta -> [ [:], fasta ] }

    SAMTOOLS_SORT (
        ch_bams_to_merge,
        ch_fasta_for_sort.first(),
        'bai'  // Generate .bai index
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    emit:
    bam      = SAMTOOLS_SORT.out.bam     // channel: [ val(meta), bam ]
    bai      = SAMTOOLS_SORT.out.bai     // channel: [ val(meta), bai ]
    versions = ch_versions               // channel: [ versions.yml ]
}
