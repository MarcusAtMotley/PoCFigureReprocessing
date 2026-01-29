/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DNA Alignment Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared Biscuit alignment for all DNA pipelines (P1, P2, P3).
    Align once, use for SNP calling, methylation, and CNV analysis.

    Input:  PE or SE FASTQs (WGS, WGEM, mTNA-DNA, HairyTNA-DNA)
    Output: Sorted BAM + BAI for downstream analysis

    Uses scatter-gather for large files:
    1. Split FASTQs into chunks using seqkit split2
    2. Align each chunk in parallel
    3. Merge BAMs back together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPLIT_FASTQ    } from '../modules/local/split_fastq/main'
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
        // SCATTER: Split FASTQs into chunks using seqkit
        //
        SPLIT_FASTQ (
            ch_reads,
            chunk_size
        )
        ch_versions = ch_versions.mix(SPLIT_FASTQ.out.versions.first())

        // Transpose to get one channel item per chunk
        // Input: [ meta, [chunk1_R1, chunk1_R2, chunk2_R1, chunk2_R2, ...] ]
        // Output: [ meta, [chunkN_R1, chunkN_R2] ] for each chunk
        ch_reads_chunked = SPLIT_FASTQ.out.reads
            .flatMap { meta, files ->
                // Group files by part number
                def chunks = [:]
                files.each { f ->
                    def match = f.name =~ /\.part_(\d+)\./
                    if (match) {
                        def part = match[0][1]
                        if (!chunks[part]) chunks[part] = []
                        chunks[part] << f
                    }
                }
                // Emit each chunk as separate item
                chunks.collect { part, chunk_files ->
                    // Sort to ensure R1 before R2
                    def sorted = chunk_files.sort { it.name }
                    def new_meta = meta.clone()
                    new_meta.chunk = part
                    [ new_meta, sorted ]
                }
            }

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
                def key = meta.id
                [ key, meta, bam ]
            }
            .groupTuple(by: 0)
            .map { key, metas, bams ->
                // Use first meta (without chunk info), flatten bams
                def clean_meta = metas[0].clone()
                clean_meta.remove('chunk')
                [ clean_meta, bams.flatten() ]
            }

        // Merge all chunks
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
