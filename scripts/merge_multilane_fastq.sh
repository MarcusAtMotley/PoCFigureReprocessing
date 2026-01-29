#!/bin/bash
# Download and merge multi-lane mTNA samples from run_motley26
# Processes ONE SAMPLE AT A TIME to minimize disk usage (~165 GiB peak)
# Run this on a VM with AWS CLI configured

set -e

S3_SRC="s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley26/results/rna_deconvolution/cutadapt"
S3_DEST="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq"

WORKDIR="${1:-fastq_merge}"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

# Process each sample one at a time to save disk space
declare -A SAMPLES=(
    ["CoB_08L_3A2_DNA-EM_S1"]="S1"
    ["CoM_08M_3A2_DNA-EM_S2"]="S2"
    ["CoB_08R_3A2_TNA-mRT-EM_S3"]="S3"
    ["CoM_08S_3A2_TNA-mRT-EM_S4"]="S4"
    ["CoB_08X_3A2_TNA-RT-EM_S5"]="S5"
    ["CoM_08Y_3A2_TNA-RT-EM_S6"]="S6"
)

for sample in CoB_08L_3A2_DNA-EM_S1 CoM_08M_3A2_DNA-EM_S2 \
              CoB_08R_3A2_TNA-mRT-EM_S3 CoM_08S_3A2_TNA-mRT-EM_S4 \
              CoB_08X_3A2_TNA-RT-EM_S5 CoM_08Y_3A2_TNA-RT-EM_S6; do

    echo ""
    echo "=========================================="
    echo "Processing: $sample"
    echo "=========================================="

    # Download this sample's files
    echo "  Downloading L001 + L002..."
    aws s3 cp "$S3_SRC/${sample}_L001_unbarcoded_R1.cutadapt.fastq" .
    aws s3 cp "$S3_SRC/${sample}_L001_unbarcoded_R2.cutadapt.fastq" .
    aws s3 cp "$S3_SRC/${sample}_L002_unbarcoded_R1.cutadapt.fastq" .
    aws s3 cp "$S3_SRC/${sample}_L002_unbarcoded_R2.cutadapt.fastq" .

    # Merge
    echo "  Merging R1..."
    cat ${sample}_L001_unbarcoded_R1.cutadapt.fastq \
        ${sample}_L002_unbarcoded_R1.cutadapt.fastq \
        > ${sample}_merged_R1.fastq

    echo "  Merging R2..."
    cat ${sample}_L001_unbarcoded_R2.cutadapt.fastq \
        ${sample}_L002_unbarcoded_R2.cutadapt.fastq \
        > ${sample}_merged_R2.fastq

    # Delete source files to free space
    echo "  Cleaning up source files..."
    rm -f ${sample}_L00*_unbarcoded_R*.cutadapt.fastq

    # Upload merged files
    echo "  Uploading merged files..."
    aws s3 cp ${sample}_merged_R1.fastq "$S3_DEST/"
    aws s3 cp ${sample}_merged_R2.fastq "$S3_DEST/"

    # Delete merged files to free space for next sample
    echo "  Cleaning up merged files..."
    rm -f ${sample}_merged_R*.fastq

    echo "  Done with $sample!"
done

echo ""
echo "=========================================="
echo "ALL DONE!"
echo "Merged files at: $S3_DEST"
echo "=========================================="
