#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PoC Figure Reprocessing Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unified pipeline for reprocessing sequencing data to generate 14 PoC figures
    comparing Motley Bio's TrinitySeq (mTNA, HairyTNA) vs traditional sequencing.

    OPTIMIZED: Shared alignments - align once per sample, use for multiple analyses.

    Pipelines:
    - P0: TrimGalore (single analyte preprocessing)
    - P1: DNA SNP calling (Revelio → LoFreq)
    - P2: DNA Methylation (Biscuit Pileup)
    - P3: CNV calling (CNVpytor)
    - P4: RNA Counts (FeatureCounts)
    - P5: RNA SNP calling (Revelio → LoFreq)

    Alignment Structure:
    - DNA samples (P1, P2, P3) → Biscuit align once → shared BAMs → P1, P2, P3
    - RNA samples (P4, P5) → STAR align once → shared BAMs → P4, P5

    Usage:
        nextflow run main.nf -profile awsbatch \
            --input samplesheet.csv \
            --outdir s3://bucket/results

    Author: Motley Bio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERGE_FASTQ   } from './pipelines/modules/local/merge_fastq/main'
include { SAMTOOLS_INDEX } from './pipelines/modules/nf-core/samtools/index/main'
include { P0_TRIMGALORE } from './pipelines/subworkflows/p0_trimgalore'
include { ALIGN_DNA     } from './pipelines/subworkflows/align_dna'
include { ALIGN_RNA     } from './pipelines/subworkflows/align_rna'
include { P1_DNA_SNP    } from './pipelines/subworkflows/p1_dna_snp'
include { P2_DNA_METH   } from './pipelines/subworkflows/p2_dna_methylation'
include { P3_CNV        } from './pipelines/subworkflows/p3_cnv'
include { P4_RNA_COUNTS } from './pipelines/subworkflows/p4_rna_counts'
include { P5_RNA_SNP    } from './pipelines/subworkflows/p5_rna_snp'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.input) {
    error "ERROR: Please provide an input samplesheet with --input"
}
if (!params.outdir) {
    error "ERROR: Please provide an output directory with --outdir"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parseSamplesheet(samplesheet_path) {
    /*
     * Parse samplesheet and handle both FASTQ and BAM inputs.
     *
     * For FASTQ input: requires fastq_1 (and optionally fastq_2)
     * For BAM input: requires bam column, bai is optional (not needed for FeatureCounts)
     *
     * Returns channel of [ meta, data ] where:
     *   - meta.input_type = 'fastq' or 'bam'
     *   - data = reads list (for FASTQ) or [bam, bai] (for BAM)
     */
    Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true, strip: true)
        .filter { row -> !row.sample.startsWith('#') }
        .map { row ->
            // Detect input type: BAM if bam column exists and is not empty
            def is_bam_input = row.bam && row.bam.trim() != ''

            def meta = [
                id             : row.sample,
                pipeline       : row.pipeline,
                cell_line      : row.cell_line ?: '',
                assay_category : row.assay_category ?: '',
                input_type     : is_bam_input ? 'bam' : 'fastq'
            ]

            if (is_bam_input) {
                // BAM input - no merging needed, already aligned
                // single_end can be specified in samplesheet (important for HairyTNA which is SE)
                meta.single_end = row.single_end?.toBoolean() ?: false
                meta.needs_merge = false
                def bam = file(row.bam.trim())
                // BAI is optional (not required for FeatureCounts, but needed for P5)
                def bai = (row.bai && row.bai.trim() != '') ? file(row.bai.trim()) : null
                meta.has_bai = (bai != null)
                [ meta, [bam, bai] ]
            } else {
                // FASTQ input
                meta.single_end = row.single_end?.toBoolean() ?: (row.fastq_2 == '' || row.fastq_2 == null)

                // Handle semicolon-separated FASTQs (multi-lane samples)
                def r1_files = row.fastq_1.split(';').collect { file(it.trim()) }
                def r2_files = row.fastq_2 ? row.fastq_2.split(';').collect { file(it.trim()) } : []

                // Mark if this sample needs lane merging
                meta.needs_merge = r1_files.size() > 1

                def reads
                if (meta.single_end) {
                    reads = r1_files
                } else {
                    // Interleave R1 and R2 for proper pairing
                    reads = [r1_files, r2_files].transpose().flatten()
                }

                [ meta, reads ]
            }
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // =========================================================================
    // Parse input samplesheet
    // =========================================================================
    ch_input = parseSamplesheet(params.input)

    // =========================================================================
    // Set up reference channels
    // =========================================================================

    ch_fasta = Channel.value([
        [ id: 'GRCh38' ],
        file(params.genome_fasta, checkIfExists: true)
    ])

    ch_fai = params.genome_fai
        ? Channel.value([ [ id: 'GRCh38' ], file(params.genome_fai, checkIfExists: true) ])
        : Channel.empty()

    ch_gtf = Channel.value([
        [ id: 'GRCh38' ],
        file(params.annotation_gtf, checkIfExists: true)
    ])

    ch_biscuit_index = params.biscuit_index
        ? Channel.value([ [ id: 'GRCh38' ], file(params.biscuit_index, checkIfExists: true) ])
        : Channel.empty()

    ch_star_index = params.star_index
        ? Channel.value([ [ id: 'GRCh38' ], file(params.star_index, checkIfExists: true) ])
        : Channel.empty()

    ch_intergenic_bed = params.intergenic_bed
        ? Channel.value(file(params.intergenic_bed, checkIfExists: true))
        : Channel.empty()

    // =========================================================================
    // Split BAM vs FASTQ input
    // =========================================================================

    ch_input.branch {
        meta, data ->
            bam_input: meta.input_type == 'bam'
            fastq_input: meta.input_type == 'fastq'
    }.set { ch_input_type }

    // =========================================================================
    // Handle BAM input - route directly to P4/P5 (pre-aligned)
    // =========================================================================

    // Split BAM input into separate bam and bai channels
    ch_bam_input = ch_input_type.bam_input.map { meta, data -> [ meta, data[0] ] }  // [ meta, bam ]
    ch_bai_input = ch_input_type.bam_input
        .filter { meta, data -> meta.has_bai }  // Only samples with BAI
        .map { meta, data -> [ meta, data[1] ] }  // [ meta, bai ]

    // P4 BAM input (FeatureCounts doesn't need BAI)
    ch_p4_bam_direct = ch_bam_input
        .filter { meta, bam -> meta.pipeline in ['P4_RNA_Counts', 'P4'] }

    // P5 BAM input (needs BAI for Revelio)
    ch_p5_bam_direct = ch_bam_input
        .filter { meta, bam -> meta.pipeline in ['P5_RNA_SNP', 'P5'] }

    // For P5 samples WITH BAI in samplesheet
    ch_p5_bai_direct = ch_bai_input
        .filter { meta, bai -> meta.pipeline in ['P5_RNA_SNP', 'P5'] }

    // For P5 samples WITHOUT BAI - generate index
    ch_p5_needs_index = ch_input_type.bam_input
        .filter { meta, data -> meta.pipeline in ['P5_RNA_SNP', 'P5'] && !meta.has_bai }
        .map { meta, data -> [ meta, data[0] ] }  // [ meta, bam ]

    SAMTOOLS_INDEX ( ch_p5_needs_index )

    // Merge BAI from samplesheet with generated BAI (for direct BAM input)
    ch_p5_bai_from_input = ch_p5_bai_direct.mix(SAMTOOLS_INDEX.out.bai)

    // =========================================================================
    // Handle FASTQ input - lane merging
    // =========================================================================

    ch_input_type.fastq_input.branch {
        meta, reads ->
            needs_merge: meta.needs_merge
            ready: !meta.needs_merge
    }.set { ch_branched }

    MERGE_FASTQ ( ch_branched.needs_merge )

    ch_merged = MERGE_FASTQ.out.reads.map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.needs_merge = false
        def read_list = reads instanceof List ? reads : [reads]
        [ new_meta, read_list ]
    }

    ch_ready = ch_branched.ready.map { meta, reads ->
        def read_list = reads instanceof List ? reads : [reads]
        [ meta, read_list ]
    }

    ch_all_reads = ch_merged.mix(ch_ready)

    // =========================================================================
    // Route FASTQ samples by analyte type (DNA vs RNA)
    // Supports comma-separated pipeline values for multi-pipeline routing
    // =========================================================================

    // Helper to parse pipeline field (supports "P1", "P1,P2", etc.)
    def parsePipelines = { pipeline_str ->
        pipeline_str.split(',').collect { it.trim() }
    }

    // Helper to check if any pipeline in list matches DNA pipelines
    def isDnaPipeline = { pipelines ->
        def dna_pipelines = ['P1_DNA_SNP', 'P1', 'P2_DNA_Meth', 'P2', 'P3_CNV', 'P3']
        pipelines.any { p -> p in dna_pipelines }
    }

    // Helper to check if any pipeline in list matches RNA pipelines
    def isRnaPipeline = { pipelines ->
        def rna_pipelines = ['P4_RNA_Counts', 'P4', 'P5_RNA_SNP', 'P5']
        pipelines.any { p -> p in rna_pipelines }
    }

    // DNA samples: P1 (SNP), P2 (Methylation), P3 (CNV)
    // De-duplicate by sample ID to align each sample only ONCE
    ch_all_reads
        .filter { meta, reads -> isDnaPipeline(parsePipelines(meta.pipeline)) }
        .map { meta, reads ->
            // Store all target pipelines in meta
            def pipelines = parsePipelines(meta.pipeline)
            def new_meta = meta.clone()
            new_meta.target_pipelines = pipelines
            [ new_meta.id, new_meta, reads ]
        }
        .groupTuple(by: 0)  // Group by sample ID (de-duplicate)
        .map { id, metas, reads_list ->
            // Merge target pipelines from all entries for this sample
            def all_pipelines = metas.collect { it.target_pipelines }.flatten().unique()
            def merged_meta = metas[0].clone()
            merged_meta.target_pipelines = all_pipelines
            [ merged_meta, reads_list[0] ]  // Use first reads (they're the same)
        }
        .set { ch_dna_reads }

    // RNA samples: P4 (Counts), P5 (SNP) - FASTQ input only
    ch_all_reads
        .filter { meta, reads -> isRnaPipeline(parsePipelines(meta.pipeline)) }
        .map { meta, reads ->
            def pipelines = parsePipelines(meta.pipeline)
            def new_meta = meta.clone()
            new_meta.target_pipelines = pipelines
            [ new_meta.id, new_meta, reads ]
        }
        .groupTuple(by: 0)
        .map { id, metas, reads_list ->
            def all_pipelines = metas.collect { it.target_pipelines }.flatten().unique()
            def merged_meta = metas[0].clone()
            merged_meta.target_pipelines = all_pipelines
            [ merged_meta, reads_list[0] ]
        }
        .set { ch_rna_reads }

    // P0: TrimGalore
    ch_all_reads
        .filter { meta, reads -> meta.pipeline in ['P0_TrimGalore', 'P0'] }
        .set { ch_p0_input }

    // =========================================================================
    // Collect versions
    // =========================================================================
    ch_versions = Channel.empty()

    // =========================================================================
    // P0: TrimGalore (preprocessing)
    // =========================================================================
    P0_TRIMGALORE ( ch_p0_input )
    ch_versions = ch_versions.mix(P0_TRIMGALORE.out.versions)

    // =========================================================================
    // SHARED DNA ALIGNMENT (Biscuit) - align once for P1, P2, P3
    // =========================================================================
    if (params.biscuit_index) {
        ALIGN_DNA (
            ch_dna_reads,
            ch_fasta,
            ch_biscuit_index
        )
        ch_versions = ch_versions.mix(ALIGN_DNA.out.versions)

        // Route aligned DNA BAMs to appropriate pipelines
        // Fan out: each BAM goes to ALL pipelines in target_pipelines
        ch_dna_bam = ALIGN_DNA.out.bam
        ch_dna_bai = ALIGN_DNA.out.bai

        // Helper to check if sample targets a specific pipeline
        def targetsP1 = { meta -> meta.target_pipelines?.any { it in ['P1_DNA_SNP', 'P1'] } }
        def targetsP2 = { meta -> meta.target_pipelines?.any { it in ['P2_DNA_Meth', 'P2'] } }
        def targetsP3 = { meta -> meta.target_pipelines?.any { it in ['P3_CNV', 'P3'] } }

        // P1: DNA SNP - filter for samples targeting P1
        ch_p1_bam = ch_dna_bam.filter { meta, bam -> targetsP1(meta) }
        ch_p1_bai = ch_dna_bai.filter { meta, bai -> targetsP1(meta) }

        P1_DNA_SNP (
            ch_p1_bam,
            ch_p1_bai,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(P1_DNA_SNP.out.versions)

        // P2: DNA Methylation - filter for samples targeting P2
        ch_p2_bam = ch_dna_bam.filter { meta, bam -> targetsP2(meta) }
        ch_p2_bai = ch_dna_bai.filter { meta, bai -> targetsP2(meta) }

        P2_DNA_METH (
            ch_p2_bam,
            ch_p2_bai,
            ch_fasta,
            ch_biscuit_index
        )
        ch_versions = ch_versions.mix(P2_DNA_METH.out.versions)

        // P3: CNV - filter for samples targeting P3
        ch_p3_bam = ch_dna_bam.filter { meta, bam -> targetsP3(meta) }
        ch_p3_bai = ch_dna_bai.filter { meta, bai -> targetsP3(meta) }

        P3_CNV (
            ch_p3_bam,
            ch_p3_bai,
            ch_fasta,
            ch_fai,
            ch_intergenic_bed
        )
        ch_versions = ch_versions.mix(P3_CNV.out.versions)
    }

    // =========================================================================
    // SHARED RNA ALIGNMENT (STAR) - align once for P4, P5 (FASTQ input only)
    // =========================================================================
    if (params.star_index) {
        ALIGN_RNA (
            ch_rna_reads,
            ch_star_index,
            ch_gtf
        )
        ch_versions = ch_versions.mix(ALIGN_RNA.out.versions)

        // Route aligned RNA BAMs to appropriate pipelines
        // Fan out: each BAM goes to ALL pipelines in target_pipelines
        ch_rna_bam = ALIGN_RNA.out.bam
        ch_rna_bai = ALIGN_RNA.out.bai

        // Helper to check if sample targets a specific pipeline
        def targetsP4 = { meta -> meta.target_pipelines?.any { it in ['P4_RNA_Counts', 'P4'] } }
        def targetsP5 = { meta -> meta.target_pipelines?.any { it in ['P5_RNA_SNP', 'P5'] } }

        // P4: RNA Counts - filter for samples targeting P4
        ch_p4_bam_aligned = ch_rna_bam.filter { meta, bam -> targetsP4(meta) }

        // P5: RNA SNP - filter for samples targeting P5
        ch_p5_bam_aligned = ch_rna_bam.filter { meta, bam -> targetsP5(meta) }
        ch_p5_bai_aligned = ch_rna_bai.filter { meta, bai -> targetsP5(meta) }

        // Mix aligned BAMs with direct BAM input
        ch_p4_bam_all = ch_p4_bam_aligned.mix(ch_p4_bam_direct)
        ch_p5_bam_all = ch_p5_bam_aligned.mix(ch_p5_bam_direct)
        ch_p5_bai_all = ch_p5_bai_aligned.mix(ch_p5_bai_from_input)

        P4_RNA_COUNTS (
            ch_p4_bam_all,
            ch_gtf
        )
        ch_versions = ch_versions.mix(P4_RNA_COUNTS.out.versions)

        P5_RNA_SNP (
            ch_p5_bam_all,
            ch_p5_bai_all,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(P5_RNA_SNP.out.versions)
    }

    // =========================================================================
    // BAM-only input (when star_index not provided)
    // =========================================================================
    if (!params.star_index) {
        // P4 with direct BAM input only
        P4_RNA_COUNTS (
            ch_p4_bam_direct,
            ch_gtf
        )
        ch_versions = ch_versions.mix(P4_RNA_COUNTS.out.versions)

        // P5 with direct BAM input only (with generated BAI if needed)
        P5_RNA_SNP (
            ch_p5_bam_direct,
            ch_p5_bai_from_input,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(P5_RNA_SNP.out.versions)
    }

    // =========================================================================
    // Collect and publish versions
    // =========================================================================
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW INTROSPECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    println ""
    println "Pipeline completed at: ${workflow.complete}"
    println "Execution status: ${workflow.success ? 'OK' : 'FAILED'}"
    println "Execution duration: ${workflow.duration}"
    println ""
}

workflow.onError {
    println "Pipeline execution stopped with error: ${workflow.errorMessage}"
}
