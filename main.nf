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
    Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true, strip: true)
        .filter { row -> !row.sample.startsWith('#') }
        .map { row ->
            def meta = [
                id             : row.sample,
                single_end     : row.single_end?.toBoolean() ?: (row.fastq_2 == '' || row.fastq_2 == null),
                pipeline       : row.pipeline,
                cell_line      : row.cell_line ?: '',
                assay_category : row.assay_category ?: ''
            ]

            def r1_files = row.fastq_1.split(';').collect { file(it.trim()) }
            def r2_files = row.fastq_2 ? row.fastq_2.split(';').collect { file(it.trim()) } : []
            meta.needs_merge = r1_files.size() > 1

            def reads
            if (meta.single_end) {
                reads = r1_files
            } else {
                reads = [r1_files, r2_files].transpose().flatten()
            }

            [ meta, reads ]
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
    // Handle lane merging
    // =========================================================================

    ch_input.branch {
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
    // Route by analyte type (DNA vs RNA)
    // =========================================================================

    // DNA samples: P1 (SNP), P2 (Methylation), P3 (CNV)
    ch_all_reads
        .filter { meta, reads ->
            meta.pipeline in ['P1_DNA_SNP', 'P1', 'P2_DNA_Meth', 'P2', 'P3_CNV', 'P3']
        }
        .set { ch_dna_reads }

    // RNA samples: P4 (Counts), P5 (SNP)
    ch_all_reads
        .filter { meta, reads ->
            meta.pipeline in ['P4_RNA_Counts', 'P4', 'P5_RNA_SNP', 'P5']
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
        ch_dna_bam = ALIGN_DNA.out.bam
        ch_dna_bai = ALIGN_DNA.out.bai

        // P1: DNA SNP - filter for P1 samples
        ch_p1_bam = ch_dna_bam.filter { meta, bam -> meta.pipeline in ['P1_DNA_SNP', 'P1'] }
        ch_p1_bai = ch_dna_bai.filter { meta, bai -> meta.pipeline in ['P1_DNA_SNP', 'P1'] }

        P1_DNA_SNP (
            ch_p1_bam,
            ch_p1_bai,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(P1_DNA_SNP.out.versions)

        // P2: DNA Methylation - filter for P2 samples
        ch_p2_bam = ch_dna_bam.filter { meta, bam -> meta.pipeline in ['P2_DNA_Meth', 'P2'] }
        ch_p2_bai = ch_dna_bai.filter { meta, bai -> meta.pipeline in ['P2_DNA_Meth', 'P2'] }

        P2_DNA_METH (
            ch_p2_bam,
            ch_p2_bai,
            ch_fasta,
            ch_biscuit_index
        )
        ch_versions = ch_versions.mix(P2_DNA_METH.out.versions)

        // P3: CNV - filter for P3 samples
        ch_p3_bam = ch_dna_bam.filter { meta, bam -> meta.pipeline in ['P3_CNV', 'P3'] }
        ch_p3_bai = ch_dna_bai.filter { meta, bai -> meta.pipeline in ['P3_CNV', 'P3'] }

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
    // SHARED RNA ALIGNMENT (STAR) - align once for P4, P5
    // =========================================================================
    if (params.star_index) {
        ALIGN_RNA (
            ch_rna_reads,
            ch_star_index,
            ch_gtf
        )
        ch_versions = ch_versions.mix(ALIGN_RNA.out.versions)

        // Route aligned RNA BAMs to appropriate pipelines
        ch_rna_bam = ALIGN_RNA.out.bam
        ch_rna_bai = ALIGN_RNA.out.bai

        // P4: RNA Counts - filter for P4 samples
        ch_p4_bam = ch_rna_bam.filter { meta, bam -> meta.pipeline in ['P4_RNA_Counts', 'P4'] }
        ch_p4_bai = ch_rna_bai.filter { meta, bai -> meta.pipeline in ['P4_RNA_Counts', 'P4'] }

        P4_RNA_COUNTS (
            ch_p4_bam,
            ch_p4_bai,
            ch_gtf
        )
        ch_versions = ch_versions.mix(P4_RNA_COUNTS.out.versions)

        // P5: RNA SNP - filter for P5 samples
        ch_p5_bam = ch_rna_bam.filter { meta, bam -> meta.pipeline in ['P5_RNA_SNP', 'P5'] }
        ch_p5_bai = ch_rna_bai.filter { meta, bai -> meta.pipeline in ['P5_RNA_SNP', 'P5'] }

        P5_RNA_SNP (
            ch_p5_bam,
            ch_p5_bai,
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
