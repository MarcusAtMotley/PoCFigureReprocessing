#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PoC Figure Reprocessing Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Unified pipeline for reprocessing sequencing data to generate 14 PoC figures
    comparing Motley Bio's TrinitySeq (mTNA, HairyTNA) vs traditional sequencing.

    Pipelines:
    - P0: TrimGalore (single analyte preprocessing)
    - P1: DNA SNP calling (Biscuit → LoFreq)
    - P2: DNA Methylation (Biscuit → Biscuit Pileup)
    - P3: CNV calling (Biscuit → CNVpytor)
    - P4: RNA Counts (STAR → FeatureCounts)
    - P5: RNA SNP calling (STAR → LoFreq)

    Usage:
        nextflow run main.nf -profile docker \
            --input samplesheet.csv \
            --outdir results

    Author: Motley Bio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERGE_FASTQ   } from './modules/local/merge_fastq/main'
include { P0_TRIMGALORE } from './subworkflows/p0_trimgalore'
include { P1_DNA_SNP    } from './subworkflows/p1_dna_snp'
include { P2_DNA_METH   } from './subworkflows/p2_dna_methylation'
include { P3_CNV        } from './subworkflows/p3_cnv'
include { P4_RNA_COUNTS } from './subworkflows/p4_rna_counts'
include { P5_RNA_SNP    } from './subworkflows/p5_rna_snp'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check required parameters
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
     * For BAM input: requires bam and bai columns
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
                meta.single_end = false  // Not applicable for BAM
                meta.needs_merge = false
                def bam = file(row.bam.trim())
                def bai = file(row.bai.trim())
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

    // Genome FASTA
    ch_fasta = Channel.value([
        [ id: 'GRCh38' ],
        file(params.genome_fasta, checkIfExists: true)
    ])

    // Genome FAI
    ch_fai = params.genome_fai
        ? Channel.value([ [ id: 'GRCh38' ], file(params.genome_fai, checkIfExists: true) ])
        : Channel.empty()

    // GTF annotation
    ch_gtf = Channel.value([
        [ id: 'GRCh38' ],
        file(params.annotation_gtf, checkIfExists: true)
    ])

    // Biscuit index (for DNA pipelines)
    ch_biscuit_index = params.biscuit_index
        ? Channel.value([ [ id: 'GRCh38' ], file(params.biscuit_index, checkIfExists: true) ])
        : Channel.empty()

    // STAR index (for RNA pipelines)
    ch_star_index = params.star_index
        ? Channel.value([ [ id: 'GRCh38' ], file(params.star_index, checkIfExists: true) ])
        : Channel.empty()

    // Intergenic BED for CNV filtering (optional)
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
    // Handle BAM input - route directly to P4/P5
    // =========================================================================

    // Split BAM input into separate bam and bai channels
    ch_bam_input = ch_input_type.bam_input.map { meta, data -> [ meta, data[0] ] }  // [ meta, bam ]
    ch_bai_input = ch_input_type.bam_input.map { meta, data -> [ meta, data[1] ] }  // [ meta, bai ]

    // =========================================================================
    // Handle FASTQ input - lane merging
    // =========================================================================

    // Split into samples needing merge vs ready to process
    ch_input_type.fastq_input.branch {
        meta, reads ->
            needs_merge: meta.needs_merge
            ready: !meta.needs_merge
    }.set { ch_branched }

    // Merge multi-lane FASTQs
    MERGE_FASTQ ( ch_branched.needs_merge )

    // Combine merged and non-merged samples
    ch_merged = MERGE_FASTQ.out.reads.map { meta, reads ->
        // Update meta to indicate merge is complete
        def new_meta = meta.clone()
        new_meta.needs_merge = false
        // Ensure reads is a list
        def read_list = reads instanceof List ? reads : [reads]
        [ new_meta, read_list ]
    }

    ch_ready = ch_branched.ready.map { meta, reads ->
        // Ensure reads is a list
        def read_list = reads instanceof List ? reads : [reads]
        [ meta, read_list ]
    }

    ch_all_reads = ch_merged.mix(ch_ready)

    // =========================================================================
    // Route FASTQ samples to pipelines
    // =========================================================================

    // P0: TrimGalore - samples needing preprocessing
    ch_all_reads
        .filter { meta, reads -> meta.pipeline == 'P0_TrimGalore' || meta.pipeline == 'P0' }
        .set { ch_p0_input }

    // P1: DNA SNP
    ch_all_reads
        .filter { meta, reads -> meta.pipeline == 'P1_DNA_SNP' || meta.pipeline == 'P1' }
        .set { ch_p1_input }

    // P2: DNA Methylation
    ch_all_reads
        .filter { meta, reads -> meta.pipeline == 'P2_DNA_Meth' || meta.pipeline == 'P2' }
        .set { ch_p2_input }

    // P3: CNV
    ch_all_reads
        .filter { meta, reads -> meta.pipeline == 'P3_CNV' || meta.pipeline == 'P3' }
        .set { ch_p3_input }

    // P4: RNA Counts - BAM input (pre-aligned samples)
    ch_p4_bam = ch_bam_input
        .filter { meta, bam -> meta.pipeline == 'P4_RNA_Counts' || meta.pipeline == 'P4' }
    ch_p4_bai = ch_bai_input
        .filter { meta, bai -> meta.pipeline == 'P4_RNA_Counts' || meta.pipeline == 'P4' }

    // P5: RNA SNP - BAM input (pre-aligned samples)
    ch_p5_bam = ch_bam_input
        .filter { meta, bam -> meta.pipeline == 'P5_RNA_SNP' || meta.pipeline == 'P5' }
    ch_p5_bai = ch_bai_input
        .filter { meta, bai -> meta.pipeline == 'P5_RNA_SNP' || meta.pipeline == 'P5' }

    // Note: FASTQ input for P4/P5 would require ALIGN_RNA first (not yet wired)
    // ch_all_reads filtered for P4/P5 is currently unused

    // =========================================================================
    // Run pipelines
    // =========================================================================

    // Collect all versions
    ch_versions = Channel.empty()

    //
    // P0: TrimGalore
    //
    P0_TRIMGALORE ( ch_p0_input )
    ch_versions = ch_versions.mix(P0_TRIMGALORE.out.versions)

    //
    // P1: DNA SNP Calling
    //
    if (params.biscuit_index) {
        P1_DNA_SNP (
            ch_p1_input,
            ch_fasta,
            ch_fai,
            ch_biscuit_index
        )
        ch_versions = ch_versions.mix(P1_DNA_SNP.out.versions)
    }

    //
    // P2: DNA Methylation
    //
    if (params.biscuit_index) {
        P2_DNA_METH (
            ch_p2_input,
            ch_fasta,
            ch_biscuit_index
        )
        ch_versions = ch_versions.mix(P2_DNA_METH.out.versions)
    }

    //
    // P3: CNV Calling
    //
    if (params.biscuit_index) {
        P3_CNV (
            ch_p3_input,
            ch_fasta,
            ch_fai,
            ch_biscuit_index,
            ch_intergenic_bed
        )
        ch_versions = ch_versions.mix(P3_CNV.out.versions)
    }

    //
    // P4: RNA Counts (BAM input)
    //
    P4_RNA_COUNTS (
        ch_p4_bam,
        ch_p4_bai,
        ch_gtf
    )
    ch_versions = ch_versions.mix(P4_RNA_COUNTS.out.versions)

    //
    // P5: RNA SNP Calling (BAM input)
    //
    P5_RNA_SNP (
        ch_p5_bam,
        ch_p5_bai,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(P5_RNA_SNP.out.versions)

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
