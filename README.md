# PoC Figure Data Reprocessing

## Goal

Reprocess sequencing data from three different library types to generate **14 PoC figures** that demonstrate the comparability and advantages of Motley Bio's TrinitySeq multi-omic assays (mTNA, HairyTNA) versus traditional single-analyte sequencing (WGS, WGEM, RNA-Seq).

---

## What is TrinitySeq?

**TrinitySeq** is Motley Bio's multi-omic sequencing technology that captures **DNA, RNA, and methylation** from a single sample in one library prep. This eliminates the need for three separate assays.

| Traditional Approach | TrinitySeq Approach |
|---------------------|---------------------|
| WGS → DNA variants + CNV | **mTNA or HairyTNA** → DNA variants + CNV |
| WGEM → Methylation | **mTNA or HairyTNA** → Methylation |
| RNA-Seq → Gene expression | **mTNA or HairyTNA** → Gene expression |
| *3 separate library preps* | *1 library prep, same data* |

### TrinitySeq Variants

| Assay | Full Name | Methylation Method | Output | Notes |
|-------|-----------|-------------------|--------|-------|
| **mTNA** | Motley TNA | Enzymatic (EM) | Paired-end | Standard TrinitySeq |
| **HairyTNA** | HairyPin TNA | Bisulfite (BS) | Single-end after collapse | Uses stem-loop adapters + UMI collapse |

The **RNA barcode** in TrinitySeq libraries allows computational separation of DNA and RNA reads after sequencing (deconvolution).

---

## Why Reprocess?

Data from different runs was processed with different pipeline versions and parameters. To make fair comparisons for the PoC figures, we need to run **all samples through identical pipelines** so differences reflect biology, not bioinformatics.

## The Problem

We have three data sources processed at different times:
- **Single Analyte** (MEDGENOME) - Raw FASTQs, never processed through our pipelines
- **mTNA** (run_motley26) - Processed, but need to verify consistency
- **HairyTNA** (run_motley34) - Processed, but different pipeline

## The Solution

1. Use **post-split intermediate FASTQs** where available (mTNA, HairyTNA)
2. Run **TrimGalore** on single analyte raw FASTQs
3. Push all samples through **5 unified pipelines** with identical parameters

---

## Quick Start

1. **Read** `DATA_REPROCESSING_PLAN.md` for full context and terminology
2. **Review** `sample_manifest.csv` for the 35 samples and their S3 paths
3. **Use** `pipeline_ready_manifest.csv` to see which FASTQs feed into which pipelines

---

## Files

| File | Description |
|------|-------------|
| `DATA_REPROCESSING_PLAN.md` | **Start here** - Comprehensive plan with data sources, pipelines, terminology, lessons learned |
| `sample_manifest.csv` | Raw FASTQ inventory (35 samples) with S3 paths and pipeline boolean flags |
| `pipeline_ready_manifest.csv` | Pipeline-specific manifest showing post-split FASTQs ready for processing |
| `processing_pipelines.md` | Detailed pipeline tool chains |
| `planning_doc.xlsx` | Original Excel planning document with Figures and Data Types sheets |

---

## Cell Lines

| Cell Line | Single Analyte | mTNA | HairyTNA |
|-----------|----------------|------|----------|
| **CoB** | ✅ WGS, WGEM, RNA | ✅ | ❌ |
| **CoM** | ✅ WGS, WGEM, RNA | ✅ | ❌ |
| **HT29** | ✅ WGS, WGEM, RNA | ❌ | ✅ (bisulfite) |

---

## Data Locations (S3)

| Source | S3 Path | Status |
|--------|---------|--------|
| Single Analyte | `s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/` | Raw - needs TrimGalore |
| mTNA (mot26) | `s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley26/results/rna_deconvolution/cutadapt/` | Post-split - ready |
| HairyTNA (mot34) | `s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/` | Post-split - ready |

**Note:** Requires AWS credentials with access to `motleybio` and `motleybio-medgenome` buckets.

---

## Pipelines

| Pipeline | Output | Tools | Inputs |
|----------|--------|-------|--------|
| **P1: DNA SNP** | VCF | Biscuit (hisat2) → Revelio → LoFreq | WGS, mTNA-DNA, HairyTNA-DNA |
| **P2: DNA Methylation** | VCF | Biscuit → Biscuit Pileup | WGEM, mTNA-DNA, HairyTNA-DNA |
| **P3: CNV** | CSV | Hisat2 → CNAnator → Intergenic Filter | WGS, mTNA-DNA, HairyTNA-DNA |
| **P4: RNA Counts** | CSV | STAR → FeatureCounts | RNA-Seq, mTNA-RNA, HairyTNA-RNA |
| **P5: RNA SNP** | VCF | STAR → Revelio → LoFreq | RNA-Seq, mTNA-RNA, HairyTNA-RNA |

**Important:** Pipelines must handle both **paired-end** (single analyte, mTNA) and **single-end** (HairyTNA) inputs!

---

## Immediate TODOs

- [ ] Run TrimGalore on 9 single analyte samples
- [ ] Build/configure 5 Nextflow mini-pipelines (PE + SE compatible)
- [ ] Create reference BED masks (intergenic, gene body, promoter, CpG islands)
- [ ] Process all samples through pipelines
- [ ] Generate PoC figures from outputs

---

## Key Terminology

| Term | Meaning |
|------|---------|
| **mTNA** | Motley TNA - enzymatic methylation, paired-end |
| **HairyTNA** | HairyPin TNA - stem-loop adapters, single-end after collapse |
| **bsTNA** | Bisulfite-converted (use these for HairyTNA - more successful) |
| **barcoded** | RNA channel (has RNA barcode) |
| **unbarcoded** | DNA channel (no RNA barcode) |

---

## PoC Figures (14 total)

| ID | Description | Data Needed |
|----|-------------|-------------|
| ID001 | DNA VAF comparison | DNA VCF |
| ID002 | DNA CNV comparison | CNV segments |
| ID003 | RNA DE comparison | Gene counts |
| ID004 | Tissue-type genes | Gene counts |
| ID005 | Exonic vs Intergenic | RSeQC |
| ID006 | Methylation by region | Meth VCF + BED masks |
| ID007 | B-cell hypomethylation | Meth VCF + B-cell BED |
| ID008-009 | Gene expression vs CNV | Counts + CNV |
| ID010 | VAF DNA vs RNA | Both VCFs |
| ID011-012 | Expression vs methylation | Counts + Meth VCF |
| ID014 | RNA deconvolution | Barcode stats |

See `DATA_REPROCESSING_PLAN.md` for full figure details.
