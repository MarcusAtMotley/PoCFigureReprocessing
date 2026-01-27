# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains planning documents and sample manifests for reprocessing sequencing data to generate 14 PoC figures comparing Motley Bio's TrinitySeq multi-omic assays (mTNA, HairyTNA) against traditional single-analyte sequencing (WGS, WGEM, RNA-Seq).

**TrinitySeq** captures DNA, RNA, and methylation from a single sample in one library prep. The goal is to reprocess all samples through identical pipelines so differences reflect biology, not bioinformatics.

## Key Concepts

### Library Types
- **mTNA**: Motley TNA - enzymatic methylation, paired-end reads
- **HairyTNA**: HairyPin TNA - bisulfite conversion, **single-end** after stem-loop collapse
- **Single Analyte**: Traditional WGS, WGEM (methylation), RNA-Seq from MEDGENOME

### Channel Terminology
- **barcoded** = RNA channel (has RNA barcode)
- **unbarcoded** = DNA channel (no RNA barcode)

### Data Sources
| Source | S3 Bucket | Status |
|--------|-----------|--------|
| Single Analyte | `s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/` | Raw - needs TrimGalore |
| mTNA (mot26) | `s3://motleybio/.../run_motley26/results/rna_deconvolution/cutadapt/` | Post-split - ready |
| HairyTNA (mot34) | `s3://motleybio/.../run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/` | Post-split - ready |

## Processing Pipelines

All pipelines must handle both paired-end (single analyte, mTNA) and single-end (HairyTNA) inputs.

| Pipeline | Output | Tools | Applies To |
|----------|--------|-------|------------|
| P1: DNA SNP | VCF | Biscuit (hisat2) → Revelio → LoFreq | WGS, mTNA-DNA, HairyTNA-DNA |
| P2: DNA Meth | VCF | Biscuit → Biscuit Pileup | WGEM, mTNA-DNA, HairyTNA-DNA |
| P3: CNV | CSV | Hisat2 → CNAnator → Intergenic Filter | WGS, mTNA-DNA, HairyTNA-DNA |
| P4: RNA Counts | CSV | STAR → FeatureCounts | RNA-Seq, mTNA-RNA, HairyTNA-RNA |
| P5: RNA SNP | VCF | STAR → Revelio → LoFreq | RNA-Seq, mTNA-RNA, HairyTNA-RNA |

## Key Files

- `sample_manifest.csv`: Raw FASTQ inventory (35 samples) with S3 paths and pipeline boolean flags
- `pipeline_ready_manifest.csv`: Pipeline-specific manifest with post-split FASTQs and processing status
- `DATA_REPROCESSING_PLAN.md`: Comprehensive plan with data sources, terminology, lessons learned

## Critical Notes

1. **HairyTNA is single-end** after stem-loop collapse - all pipelines must handle SE inputs
2. **Use bisulfite samples for HairyTNA** (bsTNA/bsRNA/bsDNA) - enzymatic conversion was less successful
3. **Revelio applied to all samples** for pipeline consistency, even non-bisulfite samples
4. **Post-split FASTQs exist** for mTNA/HairyTNA - start from intermediates, don't re-run collapse/split
5. **Cell lines**: CoB, CoM (have mTNA data), HT29 (has HairyTNA data)
