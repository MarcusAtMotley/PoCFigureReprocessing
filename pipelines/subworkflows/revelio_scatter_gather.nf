/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Revelio Scatter-Gather Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs Revelio with scatter-gather for large BAM files.

    1. Split BAM into chunks by read count
    2. Run Revelio on each chunk in parallel
    3. Merge output BAMs back together

    This avoids 16+ hour timeouts on large samples by parallelizing across nodes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPLIT_BAM       } from '../modules/local/split_bam/main'
include { REVELIO         } from '../modules/local/revelio/main'
include { SAMTOOLS_MERGE  } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX  } from '../modules/nf-core/samtools/index/main'

workflow REVELIO_SCATTER_GATHER {

    take:
    ch_bam_bai    // channel: [ val(meta), bam, bai ]
    ch_fasta      // channel: [ val(meta), fasta ]
    ch_fai        // channel: [ val(meta), fai ]

    main:
    ch_versions = Channel.empty()

    // Get chunk size from params (default 20M reads per chunk)
    // Set to 0 to disable scatter-gather
    def chunk_size = params.revelio_chunk_size ?: 20000000

    if (chunk_size > 0) {
        //
        // SCATTER: Split BAM into chunks
        //
        SPLIT_BAM (
            ch_bam_bai,
            chunk_size
        )
        ch_versions = ch_versions.mix(SPLIT_BAM.out.versions.first())

        // Flatten to get one channel item per chunk
        // Input: [ meta, [chunk1.bam, chunk2.bam, ...], [chunk1.bai, chunk2.bai, ...] ]
        // Output: [ meta, chunkN.bam, chunkN.bai ] for each chunk
        ch_bams_chunked = SPLIT_BAM.out.bams
            .flatMap { meta, bams, bais ->
                def bam_list = bams instanceof List ? bams : [bams]
                def bai_list = bais instanceof List ? bais : [bais]

                // Sort to ensure matching order
                bam_list = bam_list.sort { it.name }
                bai_list = bai_list.sort { it.name }

                // Emit each chunk as separate item with chunk number in meta
                bam_list.withIndex().collect { bam, idx ->
                    def match = bam.name =~ /\.part_(\d+)\./
                    def part = match ? match[0][1] : String.format("%03d", idx + 1)
                    def new_meta = meta.clone()
                    new_meta.chunk = part
                    [ new_meta, bam, bai_list[idx] ]
                }
            }

        //
        // MODULE: Revelio on each chunk
        //
        REVELIO (
            ch_bams_chunked,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(REVELIO.out.versions.first())

        //
        // GATHER: Group BAMs by sample ID and merge
        //
        ch_bams_grouped = REVELIO.out.bam
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
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

        // Index merged BAMs
        SAMTOOLS_INDEX ( SAMTOOLS_MERGE.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

        ch_final_bam = SAMTOOLS_MERGE.out.bam
        ch_final_bai = SAMTOOLS_INDEX.out.bai

    } else {
        //
        // SIMPLE MODE: No scatter-gather, direct Revelio
        //
        REVELIO (
            ch_bam_bai,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(REVELIO.out.versions.first())

        ch_final_bam = REVELIO.out.bam
        ch_final_bai = REVELIO.out.bai
    }

    emit:
    bam      = ch_final_bam   // channel: [ val(meta), bam ]
    bai      = ch_final_bai   // channel: [ val(meta), bai ]
    versions = ch_versions    // channel: [ versions.yml ]
}
