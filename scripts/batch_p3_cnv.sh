#!/bin/bash
# Self-contained P3 CNV pipeline for AWS Batch
# Runs CNVpytor on a markdup BAM to call copy number variants.
#
# If the input is a sorted BAM (not markdup), it will run markdup first.
#
# Usage: batch_p3_cnv.sh <SAMPLE> <BAM_PATH> [BAM_TYPE]
#   SAMPLE:    sample name
#   BAM_PATH:  S3 path to BAM file
#   BAM_TYPE:  "markdup" (default) or "sorted" (will run markdup first)
#
# Example (markdup BAM from P2):
#   batch_p3_cnv.sh CoB_08L_3A2_DNA-EM s3://bucket/p2_dna_meth/CoB_08L/CoB_08L.markdup.bam markdup
#
# Example (sorted BAM from P1, needs markdup):
#   batch_p3_cnv.sh CoB_02M_1C3_1DNA s3://bucket/p1_dna_snp/CoB_02M/CoB_02M.sorted.bam sorted
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE BAM_PATH [BAM_TYPE]}"
BAM_PATH="${2:?Missing BAM_PATH}"
BAM_TYPE="${3:-markdup}"

THREADS=$(nproc)
WORKDIR=/scratch
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
S3_REFS="s3://motleybio/Resources"
REF="$WORKDIR/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa"

echo "=========================================="
echo "AWS Batch P3 CNV Pipeline"
echo "=========================================="
echo "Sample: $SAMPLE"
echo "BAM path: $BAM_PATH"
echo "BAM type: $BAM_TYPE"
echo "Threads: $THREADS"
echo "Started: $(date)"
echo ""

# ---- Step 1: Download references ----
echo "[1/4] Downloading references..."
mkdir -p "$WORKDIR/refs" "$WORKDIR/bams" "$WORKDIR/results"

if [ ! -f "$REF" ]; then
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" "$REF"
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" "${REF}.fai"
fi
echo "References ready."

# ---- Step 2: Get analysis-ready BAM ----
ANALYSIS_BAM="$WORKDIR/bams/${SAMPLE}.analysis.bam"

echo "[2/4] Downloading BAM..."
DOWNLOADED_BAM="$WORKDIR/bams/${SAMPLE}.downloaded.bam"
aws s3 cp --quiet "$BAM_PATH" "$DOWNLOADED_BAM"

# Download or create BAI
BAI_PATH="${BAM_PATH}.bai"
if aws s3 ls "$BAI_PATH" > /dev/null 2>&1; then
    aws s3 cp --quiet "$BAI_PATH" "${DOWNLOADED_BAM}.bai"
else
    echo "  No BAI found, indexing..."
    samtools index -@ $THREADS "$DOWNLOADED_BAM"
fi

if [ "$BAM_TYPE" = "sorted" ]; then
    echo "  Running markdup on sorted BAM..."
    samtools collate -@ $THREADS -o "$WORKDIR/bams/collate.bam" "$DOWNLOADED_BAM"
    samtools fixmate -@ $THREADS -m "$WORKDIR/bams/collate.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/collate.bam"
    samtools sort -@ $THREADS -o "$WORKDIR/bams/fixmate.sorted.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/fixmate.bam"
    samtools markdup -@ $THREADS -s "$WORKDIR/bams/fixmate.sorted.bam" "$ANALYSIS_BAM" 2> "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt"
    rm -f "$WORKDIR/bams/fixmate.sorted.bam"
    samtools index -@ $THREADS "$ANALYSIS_BAM"
    rm -f "$DOWNLOADED_BAM" "${DOWNLOADED_BAM}.bai"
    echo "  Markdup done."
else
    mv "$DOWNLOADED_BAM" "$ANALYSIS_BAM"
    mv "${DOWNLOADED_BAM}.bai" "${ANALYSIS_BAM}.bai" 2>/dev/null || samtools index -@ $THREADS "$ANALYSIS_BAM"
    echo "  Using markdup BAM directly."
fi

echo "Analysis BAM ready: $(ls -lh $ANALYSIS_BAM)"

# ---- Step 3: CNVpytor ----
PYTOR_FILE="$WORKDIR/results/${SAMPLE}.pytor"
CNV_TSV="$WORKDIR/results/${SAMPLE}_cnv.tsv"

echo "[3/4] Running CNVpytor..."

# Step 3a: Read depth extraction
echo "  Extracting read depth..."
cnvpytor -root "$PYTOR_FILE" -rd "$ANALYSIS_BAM"

# Step 3b: Calculate histograms at multiple bin sizes
echo "  Calculating histograms..."
cnvpytor -root "$PYTOR_FILE" -his 1000 10000 100000

# Step 3c: Partition
echo "  Partitioning..."
cnvpytor -root "$PYTOR_FILE" -partition 1000 10000 100000

# Step 3d: Call CNVs
echo "  Calling CNVs..."
cnvpytor -root "$PYTOR_FILE" -call 1000 > "${CNV_TSV}.1000"
cnvpytor -root "$PYTOR_FILE" -call 10000 > "${CNV_TSV}.10000"
cnvpytor -root "$PYTOR_FILE" -call 100000 > "${CNV_TSV}.100000"

# Use 10000 as primary output, keep all bin sizes
cp "${CNV_TSV}.10000" "$CNV_TSV"

# Count CNV calls
CNV_COUNT=$(wc -l < "$CNV_TSV" || echo "0")
echo "CNVpytor done: $CNV_COUNT CNV calls (10kb bins)"

# Clean up BAM
rm -f "$ANALYSIS_BAM" "${ANALYSIS_BAM}.bai"

# ---- Step 4: Upload results ----
echo "[4/4] Uploading results to S3..."
aws s3 cp --quiet "$CNV_TSV" "$S3_OUTPUT/p3_cnv/${SAMPLE}/${SAMPLE}_cnv.tsv"
aws s3 cp --quiet "${CNV_TSV}.1000" "$S3_OUTPUT/p3_cnv/${SAMPLE}/${SAMPLE}_cnv_1kb.tsv"
aws s3 cp --quiet "${CNV_TSV}.10000" "$S3_OUTPUT/p3_cnv/${SAMPLE}/${SAMPLE}_cnv_10kb.tsv"
aws s3 cp --quiet "${CNV_TSV}.100000" "$S3_OUTPUT/p3_cnv/${SAMPLE}/${SAMPLE}_cnv_100kb.tsv"
aws s3 cp --quiet "$PYTOR_FILE" "$S3_OUTPUT/p3_cnv/${SAMPLE}/${SAMPLE}.pytor"

if [ -f "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt" ]; then
    aws s3 cp --quiet "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt" "$S3_OUTPUT/p3_cnv/${SAMPLE}/${SAMPLE}.markdup_metrics.txt"
fi

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE ($CNV_COUNT CNV calls)"
echo "=========================================="
echo "Finished: $(date)"
echo "Results: $S3_OUTPUT/p3_cnv/${SAMPLE}/"
