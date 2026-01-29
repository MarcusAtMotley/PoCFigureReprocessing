#!/bin/bash
# Merge already-downloaded multi-lane FASTQ files
# Use this if you already downloaded files and need to merge/upload/cleanup

set -e

S3_DEST="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq"

cd "${1:-fastq_merge}"

for sample in CoB_08L_3A2_DNA-EM_S1 CoM_08M_3A2_DNA-EM_S2 \
              CoB_08R_3A2_TNA-mRT-EM_S3 CoM_08S_3A2_TNA-mRT-EM_S4 \
              CoB_08X_3A2_TNA-RT-EM_S5 CoM_08Y_3A2_TNA-RT-EM_S6; do

    # Check if all 4 files exist for this sample
    if [ -f "${sample}_L001_unbarcoded_R1.cutadapt.fastq" ] && \
       [ -f "${sample}_L001_unbarcoded_R2.cutadapt.fastq" ] && \
       [ -f "${sample}_L002_unbarcoded_R1.cutadapt.fastq" ] && \
       [ -f "${sample}_L002_unbarcoded_R2.cutadapt.fastq" ]; then

        echo ""
        echo "=========================================="
        echo "Processing: $sample"
        echo "=========================================="

        echo "  Merging R1..."
        cat ${sample}_L001_unbarcoded_R1.cutadapt.fastq \
            ${sample}_L002_unbarcoded_R1.cutadapt.fastq > ${sample}_merged_R1.fastq

        echo "  Merging R2..."
        cat ${sample}_L001_unbarcoded_R2.cutadapt.fastq \
            ${sample}_L002_unbarcoded_R2.cutadapt.fastq > ${sample}_merged_R2.fastq

        echo "  Uploading to S3..."
        aws s3 cp ${sample}_merged_R1.fastq "$S3_DEST/"
        aws s3 cp ${sample}_merged_R2.fastq "$S3_DEST/"

        echo "  Cleaning up to free space..."
        rm -f ${sample}_L00*_unbarcoded_R*.cutadapt.fastq ${sample}_merged_R*.fastq

        echo "  Done with $sample!"
    else
        echo "Skipping $sample (files not complete yet)"
    fi
done

echo ""
echo "=========================================="
echo "Done! Merged files at: $S3_DEST"
echo "=========================================="
