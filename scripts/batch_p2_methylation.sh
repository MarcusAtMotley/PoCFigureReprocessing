#!/bin/bash
# Self-contained P2 DNA Methylation pipeline for AWS Batch
# Runs biscuit pileup on markdup BAM to extract methylation calls.
#
# Two modes:
#   BAM mode:   Download P1 sorted BAM → markdup → biscuit pileup → upload
#   FASTQ mode: Download FASTQs → biscuit align → markdup → biscuit pileup → upload
#
# Usage: batch_p2_methylation.sh <SAMPLE> <SINGLE_END> <INPUT_MODE> <INPUT_1> [INPUT_2]
#   SAMPLE:      sample name
#   SINGLE_END:  "true" | "false"
#   INPUT_MODE:  "bam" | "fastq"
#   INPUT_1:     S3 path to sorted.bam (bam mode) or semicolon-separated R1 FASTQ paths (fastq mode)
#   INPUT_2:     (fastq mode only) semicolon-separated R2 FASTQ paths
#
# Example (BAM mode):
#   batch_p2_methylation.sh CoB_08L_3A2_DNA-EM false bam s3://bucket/p1_dna_snp/CoB_08L/CoB_08L.sorted.bam
#
# Example (FASTQ mode, WGEM):
#   batch_p2_methylation.sh CoB_01W_1A3_1DNA false fastq "s3://bucket/R1.fq.gz" "s3://bucket/R2.fq.gz"
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE SINGLE_END INPUT_MODE INPUT_1 [INPUT_2]}"
SINGLE_END="${2:?Missing SINGLE_END (true/false)}"
INPUT_MODE="${3:?Missing INPUT_MODE (bam/fastq)}"
INPUT_1="${4:?Missing INPUT_1}"
INPUT_2="${5:-}"

THREADS=$(nproc)
WORKDIR=/scratch
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
S3_REFS="s3://motleybio/Resources"
REF="$WORKDIR/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BISCUIT_INDEX="$WORKDIR/refs/biscuit_index/GRCh38_full_analysis_set_plus_decoy_hla.fa"

echo "=========================================="
echo "AWS Batch P2 DNA Methylation Pipeline"
echo "=========================================="
echo "Sample: $SAMPLE"
echo "Single-end: $SINGLE_END"
echo "Input mode: $INPUT_MODE"
echo "Threads: $THREADS"
echo "Input 1: $INPUT_1"
echo "Input 2: $INPUT_2"
echo "Started: $(date)"
echo ""

# ---- Step 1: Download references ----
echo "[1/4] Downloading references..."
mkdir -p "$WORKDIR/refs/biscuit_index" "$WORKDIR/fastq" "$WORKDIR/bams" "$WORKDIR/results"

if [ ! -f "$REF" ]; then
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" "$REF"
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" "${REF}.fai"
fi

if [ ! -f "${BISCUIT_INDEX}.bis.pac" ]; then
    aws s3 cp --quiet "$S3_REFS/biscuit_reference_genome/" "$WORKDIR/refs/biscuit_index/" --recursive
fi
echo "References ready."

# ---- Step 2: Get markdup BAM ----
MARKDUP_BAM="$WORKDIR/bams/${SAMPLE}.markdup.bam"

if [ "$INPUT_MODE" = "bam" ]; then
    # BAM mode: download sorted BAM from P1, run markdup
    SORTED_BAM="$WORKDIR/bams/${SAMPLE}.sorted.bam"
    echo "[2/4] Downloading sorted BAM from S3..."
    aws s3 cp --quiet "$INPUT_1" "$SORTED_BAM"

    # Download or create BAI
    BAI_PATH="${INPUT_1}.bai"
    if aws s3 ls "$BAI_PATH" > /dev/null 2>&1; then
        aws s3 cp --quiet "$BAI_PATH" "${SORTED_BAM}.bai"
    else
        echo "  No BAI found, indexing..."
        samtools index -@ $THREADS "$SORTED_BAM"
    fi

    echo "  Running markdup..."
    samtools collate -@ $THREADS -o "$WORKDIR/bams/collate.bam" "$SORTED_BAM"
    samtools fixmate -@ $THREADS -m "$WORKDIR/bams/collate.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/collate.bam"
    samtools sort -@ $THREADS -o "$WORKDIR/bams/fixmate.sorted.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/fixmate.bam"
    samtools markdup -@ $THREADS -s "$WORKDIR/bams/fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt"
    rm -f "$WORKDIR/bams/fixmate.sorted.bam"
    samtools index -@ $THREADS "$MARKDUP_BAM"
    rm -f "$SORTED_BAM" "${SORTED_BAM}.bai"
    echo "  Markdup done."

elif [ "$INPUT_MODE" = "fastq" ]; then
    # FASTQ mode: download FASTQs, align with biscuit, then markdup
    echo "[2/4] Downloading FASTQs and aligning..."

    FQ1="$WORKDIR/fastq/${SAMPLE}_R1.fastq"
    FQ2="$WORKDIR/fastq/${SAMPLE}_R2.fastq"

    # Download R1 files
    IDX=0
    IFS=';' read -ra R1_URLS <<< "$INPUT_1"
    for url in "${R1_URLS[@]}"; do
        LOCAL="$WORKDIR/fastq/r1_${IDX}"
        echo "  Downloading R1[$IDX]: $url"
        aws s3 cp --quiet "$url" "$LOCAL"
        IDX=$((IDX+1))
    done

    # Concatenate R1 lanes
    if [ ${#R1_URLS[@]} -gt 1 ]; then
        cat "$WORKDIR/fastq"/r1_* > "$FQ1"
        rm -f "$WORKDIR/fastq"/r1_*
    else
        mv "$WORKDIR/fastq/r1_0" "$FQ1"
    fi

    # Download R2 files (if paired-end)
    if [ "$SINGLE_END" = "false" ] && [ -n "$INPUT_2" ]; then
        IDX=0
        IFS=';' read -ra R2_URLS <<< "$INPUT_2"
        for url in "${R2_URLS[@]}"; do
            LOCAL="$WORKDIR/fastq/r2_${IDX}"
            echo "  Downloading R2[$IDX]: $url"
            aws s3 cp --quiet "$url" "$LOCAL"
            IDX=$((IDX+1))
        done
        if [ ${#R2_URLS[@]} -gt 1 ]; then
            cat "$WORKDIR/fastq"/r2_* > "$FQ2"
            rm -f "$WORKDIR/fastq"/r2_*
        else
            mv "$WORKDIR/fastq/r2_0" "$FQ2"
        fi
    fi

    echo "  FASTQs ready."
    ls -lh "$WORKDIR/fastq/"

    # Biscuit alignment
    SORTED_BAM="$WORKDIR/bams/${SAMPLE}.sorted.bam"
    echo "  Running Biscuit alignment..."
    if [ "$SINGLE_END" = "true" ]; then
        biscuit align -@ $THREADS "$BISCUIT_INDEX" "$FQ1" \
            | samtools sort -@ 8 -o "$SORTED_BAM" -
    else
        biscuit align -@ $THREADS "$BISCUIT_INDEX" "$FQ1" "$FQ2" \
            | samtools sort -@ 8 -o "$SORTED_BAM" -
    fi
    samtools index -@ $THREADS "$SORTED_BAM"
    echo "  Alignment done."

    # Upload sorted BAM as checkpoint (for P3 reuse)
    echo "  Uploading alignment checkpoint..."
    aws s3 cp --quiet "$SORTED_BAM" "$S3_OUTPUT/p2_dna_meth/${SAMPLE}/${SAMPLE}.sorted.bam"
    aws s3 cp --quiet "${SORTED_BAM}.bai" "$S3_OUTPUT/p2_dna_meth/${SAMPLE}/${SAMPLE}.sorted.bam.bai"

    # Free FASTQs
    rm -f "$WORKDIR/fastq/"*

    # Markdup
    echo "  Running markdup..."
    samtools collate -@ $THREADS -o "$WORKDIR/bams/collate.bam" "$SORTED_BAM"
    samtools fixmate -@ $THREADS -m "$WORKDIR/bams/collate.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/collate.bam"
    samtools sort -@ $THREADS -o "$WORKDIR/bams/fixmate.sorted.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/fixmate.bam"
    samtools markdup -@ $THREADS -s "$WORKDIR/bams/fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt"
    rm -f "$WORKDIR/bams/fixmate.sorted.bam"
    samtools index -@ $THREADS "$MARKDUP_BAM"
    rm -f "$SORTED_BAM" "${SORTED_BAM}.bai"
    echo "  Markdup done."

else
    echo "ERROR: Unknown INPUT_MODE '$INPUT_MODE'. Expected 'bam' or 'fastq'."
    exit 1
fi

echo "Markdup BAM ready: $(ls -lh $MARKDUP_BAM)"

# ---- Step 3: Biscuit pileup (methylation extraction) ----
METH_VCF="$WORKDIR/results/${SAMPLE}.methylation.vcf"
echo "[3/4] Running biscuit pileup (methylation extraction)..."
biscuit pileup -q 20 -@ $THREADS "$BISCUIT_INDEX" "$MARKDUP_BAM" -o "$METH_VCF"

bgzip -c "$METH_VCF" > "${METH_VCF}.gz"
tabix -p vcf "${METH_VCF}.gz"

# Count methylation sites
METH_COUNT=$(grep -vc '^#' "$METH_VCF" || echo "0")
echo "Methylation extraction done: $METH_COUNT sites"
rm -f "$METH_VCF"

# ---- Step 4: Upload results ----
echo "[4/4] Uploading results to S3..."
aws s3 cp --quiet "${METH_VCF}.gz" "$S3_OUTPUT/p2_dna_meth/${SAMPLE}/${SAMPLE}.methylation.vcf.gz"
aws s3 cp --quiet "${METH_VCF}.gz.tbi" "$S3_OUTPUT/p2_dna_meth/${SAMPLE}/${SAMPLE}.methylation.vcf.gz.tbi"
aws s3 cp --quiet "$MARKDUP_BAM" "$S3_OUTPUT/p2_dna_meth/${SAMPLE}/${SAMPLE}.markdup.bam"
aws s3 cp --quiet "${MARKDUP_BAM}.bai" "$S3_OUTPUT/p2_dna_meth/${SAMPLE}/${SAMPLE}.markdup.bam.bai"
aws s3 cp --quiet "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt" "$S3_OUTPUT/p2_dna_meth/${SAMPLE}/${SAMPLE}.markdup_metrics.txt"

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE ($METH_COUNT methylation sites)"
echo "=========================================="
echo "Finished: $(date)"
echo "Results: $S3_OUTPUT/p2_dna_meth/${SAMPLE}/"
