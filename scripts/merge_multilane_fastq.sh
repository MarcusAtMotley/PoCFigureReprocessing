#!/bin/bash
# Download and merge multi-lane mTNA samples from run_motley26
# Run this on a VM with AWS CLI configured

set -e

S3_SRC="s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley26/results/rna_deconvolution/cutadapt"
S3_DEST="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq"

WORKDIR="${1:-fastq_merge}"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

echo "=== Downloading multi-lane files from $S3_SRC ==="

# Download all multi-lane files (S1-S6, L001+L002)
aws s3 cp "$S3_SRC/" . --recursive --exclude "*" \
    --include "*_S1_L00*_unbarcoded_R*.cutadapt.fastq" \
    --include "*_S2_L00*_unbarcoded_R*.cutadapt.fastq" \
    --include "*_S3_L00*_unbarcoded_R*.cutadapt.fastq" \
    --include "*_S4_L00*_unbarcoded_R*.cutadapt.fastq" \
    --include "*_S5_L00*_unbarcoded_R*.cutadapt.fastq" \
    --include "*_S6_L00*_unbarcoded_R*.cutadapt.fastq"

echo "=== Merging lanes ==="

# Merge each sample (L001 + L002 -> merged)
for sample in CoB_08L_3A2_DNA-EM_S1 CoM_08M_3A2_DNA-EM_S2 \
              CoB_08R_3A2_TNA-mRT-EM_S3 CoM_08S_3A2_TNA-mRT-EM_S4 \
              CoB_08X_3A2_TNA-RT-EM_S5 CoM_08Y_3A2_TNA-RT-EM_S6; do

    echo "  Merging $sample R1..."
    cat ${sample}_L001_unbarcoded_R1.cutadapt.fastq \
        ${sample}_L002_unbarcoded_R1.cutadapt.fastq \
        > ${sample}_merged_R1.fastq

    echo "  Merging $sample R2..."
    cat ${sample}_L001_unbarcoded_R2.cutadapt.fastq \
        ${sample}_L002_unbarcoded_R2.cutadapt.fastq \
        > ${sample}_merged_R2.fastq
done

echo "=== Uploading merged files to $S3_DEST ==="

aws s3 cp . "$S3_DEST/" --recursive --exclude "*" --include "*_merged_R*.fastq"

echo "=== Done! ==="
echo "Merged files uploaded to: $S3_DEST"
echo ""
echo "Files created:"
ls -lh *_merged_R*.fastq
