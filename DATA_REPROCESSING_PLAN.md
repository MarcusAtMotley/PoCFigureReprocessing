# Motley Bio Data Reprocessing Plan

## Overview

This document captures the comprehensive plan for reprocessing sequencing data to generate comparable figures across different library types (Single Analyte, mTNA, HairyTNA) for the PoC figures.

**Date:** 2026-01-27
**Cell Lines of Interest:** CoM, CoB, HT29

---

## Data Sources

### 1. Single Analyte (MEDGENOME P2008501) - Deep Sequencing

**Location:** `s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/`

| Cell Line | WGS (DNA SNPs/CNV) | WGEM (Methylation) | RNA-Seq |
|-----------|-------------------|---------------------|---------|
| **CoB** | `CoB_02M_1C3_1DNA` | `CoB_01W_1A3_1DNA` | `CoB_02W_1A3_1RNA` |
| **CoM** | `CoM_02K_1C3_1DNA` | `CoM_01T_1A3_1DNA` | `CoM_02V_1A3_1RNA` |
| **HT29** | `HT2_02N_1B3_1DNA` | `HT2_01Z_1A3_1DNA` | `HT2_02T_1A3_1RNA` |

**Status:** Raw FASTQs only - need TrimGalore then pipeline processing
**Total Size:** ~1.3 TB

**Naming Convention:**
- `1A3` = WGEM (enzymatic methylation)
- `1C3` / `1B3` = WGS (standard DNA)
- `1RNA` = RNA-Seq

---

### 2. Deep mTNA (run_motley26)

**Location:** `s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley26/results/`

| Cell Line | Sample | Library Type | Notes |
|-----------|--------|--------------|-------|
| **CoB** | `CoB_08L_3A2_DNA-EM` | DNA control | 2 lanes |
| **CoB** | `CoB_08R_3A2_TNA-mRT-EM` | mTNA (mRT) | 2 lanes |
| **CoB** | `CoB_08X_3A2_TNA-RT-EM` | mTNA (RT) | 2 lanes |
| **CoB** | `CoB_09D_3A2_RNA-mRT-EM` | RNA (mRT) | 2 lanes |
| **CoB** | `CoB_09J_3A2_RNA-RT-EM` | RNA (RT) | 2 lanes |
| **CoM** | `CoM_08M_3A2_DNA-EM` | DNA control | 2 lanes |
| **CoM** | `CoM_08S_3A2_TNA-mRT-EM` | mTNA (mRT) | 2 lanes |
| **CoM** | `CoM_08Y_3A2_TNA-RT-EM` | mTNA (RT) | 2 lanes |
| **CoM** | `CoM_09E_3A2_RNA-mRT-EM` | RNA (mRT) | 2 lanes |
| **CoM** | `CoM_09K_3A2_RNA-RT-EM` | RNA (RT) | 2 lanes |

**Status:** Post-split FASTQs available in `rna_deconvolution/cutadapt/`
**Total Size:** ~175 GB

**Key Paths:**
- Raw FASTQs: `results/bclconvert/output/`
- **Post-split FASTQs:** `results/rna_deconvolution/cutadapt/`
  - `*_barcoded_R1.cutadapt.fastq` / `*_barcoded_R2.cutadapt.fastq` → RNA channel
  - `*_unbarcoded_R1.cutadapt.fastq` / `*_unbarcoded_R2.cutadapt.fastq` → DNA channel

---

### 3. Deep HairyTNA (run_motley34) - Bisulfite Only

**Location:** `s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/`

| Cell Line | Sample | Library Type | Notes |
|-----------|--------|--------------|-------|
| **HT29** | `HT29_21S_3A2_bsTNA-HP` | bsTNA | Replicate 1 |
| **HT29** | `HT29_21T_3A2_bsTNA-HP` | bsTNA | Replicate 2 |
| **HT29** | `HT29_21U_3A2_bsRNA-HP` | bsRNA | Replicate 1 |
| **HT29** | `HT29_21V_3A2_bsRNA-HP` | bsRNA | Replicate 2 |
| **HT29** | `HT29_21W_3A2_bsDNA-HP` | bsDNA | Replicate 1 |
| **HT29** | `HT29_21X_3A2_bsDNA-HP` | bsDNA | Replicate 2 |

**Important:** HairyTNA uses **bisulfite** (bs) conversion, which was more successful than enzymatic (em) for these samples. Ignore the emTNA/emRNA/emDNA samples.

**Status:** Post-collapse/split FASTQs available
**Total Size:** ~16 GB

**Key Paths:**
- Raw FASTQs: `bcl_and_rawfastqs/fastq/`
- Post-collapse consensus: `seqera_HP_output/stem_collapse/*_consensus.fastq.gz`
- **Post-split FASTQs:** `seqera_HP_output/rna_deconvolution/cutadapt/`
  - `*_barcoded.cutadapt.fastq` → RNA channel (SINGLE-END!)
  - `*_unbarcoded.cutadapt.fastq` → DNA channel (SINGLE-END!)

**Critical Note:** HairyTNA FASTQs are **SINGLE-END** after stem-loop collapse, not paired-end like mTNA and single analyte!

---

## Library Type Terminology

| Term | Meaning |
|------|---------|
| **mTNA** | Motley TNA - enzymatic methylation conversion |
| **HairyTNA** | HairyPin TNA - uses stem-loop adapters, requires UMI collapse |
| **bsTNA/bsRNA/bsDNA** | Bisulfite-converted samples (more successful for HairyTNA) |
| **emTNA/emRNA/emDNA** | Enzymatic-converted samples (less successful for HairyTNA) |
| **mRT vs RT** | Different reverse transcriptase variants in library prep |
| **HP vs pHP** | HairyPin vs protein-barcoded HairyPin (ignore pHP for now) |
| **barcoded** | RNA channel (has sequencing barcode for RNA identification) |
| **unbarcoded** | DNA channel (no RNA barcode) |

---

## Unified Processing Pipelines

To ensure comparability, all samples go through the same downstream pipelines after appropriate preprocessing.

### Pipeline Overview

| Pipeline | Output | Input Samples | Tools |
|----------|--------|---------------|-------|
| **P0: TrimGalore** | Trimmed FASTQ | Single analyte only | TrimGalore |
| **P1: DNA SNP** | VCF | WGS, mTNA-DNA, HairyTNA-DNA | Biscuit (hisat2) → Revelio → LoFreq |
| **P2: DNA Methylation** | VCF | WGEM, mTNA-DNA, HairyTNA-DNA | Biscuit → Biscuit Pileup |
| **P3: CNV** | CSV | WGS, mTNA-DNA, HairyTNA-DNA | Hisat2 → CNAnator → Intergenic Filter |
| **P4: RNA Counts** | CSV | RNA-Seq, mTNA-RNA, HairyTNA-RNA | STAR → FeatureCounts |
| **P5: RNA SNP** | VCF | RNA-Seq, mTNA-RNA, HairyTNA-RNA | STAR → Revelio → LoFreq |

### Pipeline Details

#### P0: TrimGalore (Single Analyte Preprocessing)
```
Raw PE FASTQ → TrimGalore → Trimmed PE FASTQ
```
- Only needed for single analyte samples
- mTNA and HairyTNA already have trimmed/processed intermediates

#### P1: DNA SNP Calling
```
[PE or SE FASTQ] → Biscuit Align (hisat2 mode) → Revelio BS reversal → LoFreq → VCF
```
- Revelio applied even to non-BS samples for consistency
- Handles both PE (single analyte, mTNA) and SE (HairyTNA) inputs

#### P2: DNA Methylation
```
[PE or SE FASTQ] → Biscuit Align → Biscuit Pileup → Methylation VCF
```

#### P3: CNV Calling
```
[PE or SE FASTQ] → Hisat2 Align → CNAnator → Intergenic Filter → CNV CSV
```
- CNAnator chosen over ichorCNA for unified approach
- Intergenic filtering removes noise from non-coding regions

#### P4: RNA Gene Counts
```
[PE or SE FASTQ] → STAR Align → FeatureCounts → Counts CSV
```

#### P5: RNA SNP Calling
```
[PE or SE FASTQ] → STAR Align → Revelio BS reversal → LoFreq → VCF
```

---

## Sample × Pipeline Matrix

### Which samples need which pipelines:

| Sample Type | P1 DNA SNP | P2 DNA Meth | P3 CNV | P4 RNA Counts | P5 RNA SNP |
|-------------|------------|-------------|--------|---------------|------------|
| **WGS** | ✅ | ❌ | ✅ | ❌ | ❌ |
| **WGEM** | ❌ | ✅ | ❌ | ❌ | ❌ |
| **RNA-Seq** | ❌ | ❌ | ❌ | ✅ | ✅ |
| **mTNA** (DNA ch) | ✅ | ✅ | ✅ | ❌ | ❌ |
| **mTNA** (RNA ch) | ❌ | ❌ | ❌ | ✅ | ✅ |
| **HairyTNA** (DNA ch) | ✅ | ✅ | ✅ | ❌ | ❌ |
| **HairyTNA** (RNA ch) | ❌ | ❌ | ❌ | ✅ | ✅ |

---

## Data Readiness Summary

### What's Ready vs Needs Processing

| Status | Count | Description |
|--------|-------|-------------|
| **ready** | 40 pipeline-runs | Post-split mTNA and HairyTNA FASTQs |
| **TrimGalore needed** | 15 pipeline-runs | Single analyte raw FASTQs (9 unique samples) |

### Read Type Considerations

| Read Type | Sources | Pipeline Implications |
|-----------|---------|----------------------|
| **Paired-End (PE)** | Single analyte, mTNA | Standard PE alignment |
| **Single-End (SE)** | HairyTNA | Need SE-compatible alignment params |

All pipelines must handle both PE and SE inputs!

---

## Figure Requirements (from Excel planning doc)

| Fig ID | Description | Data Needed |
|--------|-------------|-------------|
| ID001 | DNA VAF comparison | VCF (DNA channel) |
| ID002 | DNA CNV comparison | CNV seg files |
| ID003 | RNA DE comparison | FeatureCounts/RSEM |
| ID004 | Tissue-type genes | Gene count files |
| ID005 | Exonic vs Intergenic | RSeQC outputs |
| ID006 | Methylation regions | Meth VCF + BED masks |
| ID007 | B-cell hypomethylation | Meth VCF + B-cell BED |
| ID008 | Gene expr vs CNV (trinity) | Counts + CNV |
| ID009 | Gene expr vs CNV (single) | Counts + CNV |
| ID010 | VAF DNA vs RNA | VCF both channels |
| ID011 | Gene expr vs methylation (trinity) | Counts + Meth VCF |
| ID012 | Gene expr vs methylation (single) | Counts + Meth VCF |
| ID014 | RNA deconvolution | Custom barcode stats |

---

## Required Reference Files

1. **Genome Reference:** GRCh38
   - Biscuit index
   - Hisat2 index
   - STAR index

2. **Gene Annotation:** Gencode GTF (for FeatureCounts)

3. **BED Masks:**
   - Intergenic regions (for CNV filtering)
   - Gene bodies (for methylation analysis)
   - Promoters (TSS ± 2kb)
   - CpG islands
   - B-cell hypomethylation sites

---

## Files Created

| File | Description |
|------|-------------|
| `sample_manifest.csv` | Complete sample inventory with raw FASTQ paths and pipeline flags |
| `pipeline_ready_manifest.csv` | Pipeline-specific manifest with post-split intermediate FASTQs |
| `processing_pipelines.md` | Pipeline documentation |
| `DATA_REPROCESSING_PLAN.md` | This file |

---

## Next Steps

### Phase 1: Single Analyte Preprocessing
1. Run TrimGalore on 9 single analyte samples
2. Output trimmed FASTQs to consistent location

### Phase 2: Build Mini-Nextflow Pipelines
Create 5 small Nextflow workflows that:
- Accept both PE and SE inputs
- Use consistent tool versions
- Output to standardized locations

### Phase 3: Run All Samples Through Pipelines
- ~55 total pipeline runs across all samples
- Parallelize where possible

### Phase 4: Generate Region Masks
- Create BED files for methylation regional analysis
- Source from UCSC/Gencode annotations

---

## Key Lessons Learned

1. **HairyTNA is single-end** after stem-loop collapse - pipelines must handle both PE and SE

2. **Bisulfite > Enzymatic for HairyTNA** - use bsTNA/bsRNA/bsDNA samples, not emTNA

3. **Post-split FASTQs exist** - don't need to re-run collapse/split, can start from intermediates

4. **Channel terminology:**
   - "barcoded" = RNA channel
   - "unbarcoded" = DNA channel

5. **Revelio for consistency** - apply BS reversal even to non-BS samples for comparable variant calling

6. **TippingPoint samples are different cell types** - TP-B6, TP-CN2, TP-Par are NOT CoM/CoB/HT29

7. **Lane merging decision** - mot26 samples have L001/L002 lanes, may need to merge or process separately
