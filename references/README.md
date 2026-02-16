# Reference BED Files for Post-hoc CNV Filtering

## Purpose
TrinitySeq captures DNA and RNA in a single library prep. IchorCNA CNV calls on the DNA channel may contain artifacts from RNA bleedthrough in genic regions. These BED files enable post-hoc filtering with `bedtools intersect`.

## Files

| File | Regions | Count | Coverage |
|------|---------|-------|----------|
| `intergenic_regions.hg38.bed` | Between genes | 32,661 | 1.31 Gbp (42.3%) |
| `genic_regions.hg38.bed` | Merged gene bodies | 32,636 | 1.78 Gbp (57.7%) |

Also uploaded to: `s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/references/`

## How They Were Generated

**Source:** Ensembl GRCh38 release 112 transcript BED (`s3://motleybio/Resources/GTF/Homo_sapiens.GRCh38.112.chr_label.bed`)

**Steps:**
```bash
# 1. Extract chr/start/end from 254K transcript-level BED, merge overlapping genes
cut -f1-3 ensembl_genes.bed | sort -k1,1 -k2,2n | bedtools merge -i - > genic_regions.bed

# 2. Filter to standard chromosomes (chr1-22, X, Y, M), fix chrMT→chrM
grep -E '^chr([0-9]+|X|Y|MT)\b' genic_regions.bed | grep -v '_' | sed 's/^chrMT/chrM/' > genic_std.bed

# 3. Sort by genome order
bedtools sort -i genic_std.bed -g genome.sizes > genic_regions.hg38.bed

# 4. Complement against genome to get intergenic regions
bedtools complement -i genic_regions.hg38.bed -g genome.sizes > intergenic_regions.hg38.bed
```

**Genome sizes** from: `s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai`

## Usage

Filter CNV calls to intergenic regions only (removes RNA bleedthrough artifacts):
```bash
bedtools intersect -a cnv_calls.bed -b intergenic_regions.hg38.bed > cnv_intergenic_only.bed
```

Or annotate which CNV segments overlap genic regions:
```bash
bedtools intersect -a cnv_calls.bed -b genic_regions.hg38.bed -wa -wb > cnv_with_gene_overlap.bed
```

---
Generated: Feb 16, 2026 | bedtools v2.31.1
