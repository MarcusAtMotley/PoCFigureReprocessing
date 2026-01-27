# Unified Processing Pipelines for PoC Figures

## Overview

To ensure comparability across all library types (Single Analyte, mTNA, HairyTNA), we need to process all samples through **unified pipelines** that use the same tools and parameters.

**Note:** `*` in the original spreadsheet indicates steps that are technically unnecessary but included for consistency.

---

## Pipeline Summary

### Pipeline 1: DNA SNP Calling (VCF)
**Output:** VCF with DNA SNPs
**Applies to:** WGS, mTNA (DNA channel), HairyTNA (DNA channel)

```
FASTQ → Biscuit Align [hisat2 mode] → Revelio BS reversal → LoFreq Call → VCF
```

**Tools:**
- Biscuit (alignment with hisat2 mode)
- Revelio (bisulfite conversion reversal - for consistency even on non-BS samples)
- LoFreq (variant calling)

---

### Pipeline 2: DNA Methylation (VCF)
**Output:** VCF with methylation calls
**Applies to:** WGEM, mTNA (DNA channel), HairyTNA (DNA channel)

```
FASTQ → Filter → Biscuit Align → Biscuit Pileup → Methylation VCF
```

**Tools:**
- TrimGalore/Cutadapt (filtering)
- Biscuit (alignment + pileup)

---

### Pipeline 3: DNA Copy Number (CNV CSV)
**Output:** CNV segments with intergenic filtering
**Applies to:** WGS, mTNA (DNA channel), HairyTNA (DNA channel)

```
FASTQ → Hisat2 Align → CNAnator Call → Intergenic Filter → CNV CSV
```

**Tools:**
- Hisat2 (alignment)
- CNAnator (CNV calling)
- Custom intergenic filter (BED mask)

---

### Pipeline 4: RNA Gene Counts (Counts CSV)
**Output:** Gene count matrix
**Applies to:** RNA-Seq, mTNA (RNA channel), HairyTNA (RNA channel)

```
FASTQ → STAR Align → FeatureCounts → Counts CSV
```

**Tools:**
- STAR (alignment)
- FeatureCounts (quantification)

---

### Pipeline 5: RNA SNP Calling (VCF)
**Output:** VCF with RNA SNPs
**Applies to:** RNA-Seq, mTNA (RNA channel), HairyTNA (RNA channel)

```
FASTQ → STAR Align → Revelio BS reversal → LoFreq Call → VCF
```

**Tools:**
- STAR (alignment)
- Revelio (for consistency)
- LoFreq (variant calling)

---

## Pre-Processing Steps by Library Type

### Single Analyte (WGS, WGEM, RNA-Seq)
```
Raw FASTQ → TrimGalore → Pipeline(s)
```
No splitting required - direct to pipeline.

### mTNA
```
Raw FASTQ → TrimGalore → Barcode Split → DNA Channel → DNA Pipelines (1, 2, 3)
                                       → RNA Channel → RNA Pipelines (4, 5)
```

### HairyTNA
```
Raw FASTQ → TrimGalore → UMI Collapse → Barcode Split → DNA Channel → DNA Pipelines (1, 2, 3)
                                                       → RNA Channel → RNA Pipelines (4, 5)
```

---

## Sample × Pipeline Matrix

| Sample Type | Pipeline 1 (DNA SNP) | Pipeline 2 (DNA Meth) | Pipeline 3 (CNV) | Pipeline 4 (RNA Counts) | Pipeline 5 (RNA SNP) |
|-------------|---------------------|----------------------|------------------|------------------------|---------------------|
| **WGS** | ✅ | ❌ | ✅ | ❌ | ❌ |
| **WGEM** | ❌ | ✅ | ❌ | ❌ | ❌ |
| **RNA-Seq** | ❌ | ❌ | ❌ | ✅ | ✅ |
| **mTNA** | ✅ (DNA ch) | ✅ (DNA ch) | ✅ (DNA ch) | ✅ (RNA ch) | ✅ (RNA ch) |
| **HairyTNA** | ✅ (DNA ch) | ✅ (DNA ch) | ✅ (DNA ch) | ✅ (RNA ch) | ✅ (RNA ch) |

---

## Required Reference Files

1. **Genome Reference:** GRCh38 (with Biscuit index, Hisat2 index, STAR index)
2. **Gene Annotation:** Gencode GTF (for FeatureCounts)
3. **Intergenic Mask:** BED file for CNV filtering
4. **Gene Body Mask:** BED file for methylation region analysis
5. **Promoter Mask:** BED file (TSS ± 2kb)
6. **CpG Islands Mask:** BED file
7. **B-cell Hypomethylation Sites:** BED file

---

## Processing Priority

### Phase 1: Single Analyte (9 samples)
- Establishes "gold standard" comparison data
- ~1.3 TB of data to process

### Phase 2: mTNA Reprocessing (20 samples from mot26)
- May already have some outputs, but reprocess for consistency
- ~175 GB of data

### Phase 3: HairyTNA Reprocessing (6 samples from mot34)
- Bisulfite samples only (bsTNA, bsRNA, bsDNA)
- ~16 GB of data

---

## Notes

- Revelio BS reversal is applied even to non-bisulfite samples for pipeline consistency
- CNAnator chosen over ichorCNA for unified approach
- All samples should use identical tool versions and parameters
