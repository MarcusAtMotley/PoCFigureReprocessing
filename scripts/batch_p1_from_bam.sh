#!/bin/bash
# P1 DNA SNP pipeline — resume from sorted BAM
# Skips alignment, picks up at markdup → calmd → revelio → bcftools
#
# Use when:
#   - Alignment completed but post-alignment OOM'd (resubmit with more memory)
#   - Running WGEM through P1 using markdup BAMs from P2 (set --start-from markdup)
#
# Usage: batch_p1_from_bam.sh <SAMPLE> <BAM_S3_PATH> [BAM_TYPE]
#   SAMPLE:       sample name
#   BAM_S3_PATH:  S3 path to sorted.bam or markdup.bam
#   BAM_TYPE:     "sorted" (default) — runs markdup first
#                 "markdup" — skips markdup, starts at calmd
#
# Example (resume from sorted BAM after OOM):
#   batch_p1_from_bam.sh HT29_02N_1B3_1DNA s3://bucket/p1_dna_snp/HT29_02N/HT29_02N.sorted.bam sorted
#
# Example (WGEM from P2 markdup BAM):
#   batch_p1_from_bam.sh CoB_01W_1A3_1DNA s3://bucket/p2_dna_meth/CoB_01W/CoB_01W.markdup.bam markdup
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE BAM_S3_PATH [BAM_TYPE]}"
BAM_S3="${2:?Missing BAM_S3_PATH}"
BAM_TYPE="${3:-sorted}"

THREADS=$(nproc)
WORKDIR=/scratch
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
S3_REFS="s3://motleybio/Resources"
REF="$WORKDIR/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa"

echo "=========================================="
echo "P1 DNA SNP — Resume from BAM"
echo "=========================================="
echo "Sample: $SAMPLE"
echo "BAM source: $BAM_S3"
echo "BAM type: $BAM_TYPE"
echo "Threads: $THREADS"
echo "Started: $(date)"
echo ""

# ---- Clean scratch from prior runs ----
rm -rf "$WORKDIR/bams" "$WORKDIR/results" "$WORKDIR/revelio_tmp"
mkdir -p "$WORKDIR/refs" "$WORKDIR/bams" "$WORKDIR/results"

# ---- Download reference ----
echo "[1] Downloading reference..."
if [ ! -f "$REF" ]; then
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" "$REF"
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" "${REF}.fai"
fi
echo "Reference ready."

# ---- Download BAM ----
echo "[2] Downloading BAM..."
INPUT_BAM="$WORKDIR/bams/${SAMPLE}.input.bam"
aws s3 cp --quiet "$BAM_S3" "$INPUT_BAM"
# Download BAI if it exists
aws s3 cp --quiet "${BAM_S3}.bai" "${INPUT_BAM}.bai" 2>/dev/null || true
echo "BAM ready: $(ls -lh $INPUT_BAM)"

# Index if BAI wasn't available
if [ ! -f "${INPUT_BAM}.bai" ]; then
    echo "  Indexing BAM..."
    samtools index -@ $THREADS "$INPUT_BAM"
fi

# ---- Mark duplicates (skip if BAM_TYPE=markdup) ----
MARKDUP_BAM="$WORKDIR/bams/${SAMPLE}.markdup.bam"

if [ "$BAM_TYPE" = "sorted" ]; then
    echo "[3] Running samtools markdup..."
    echo "  collate → fixmate → sort → markdup"
    echo "  Using reduced threads for sort to limit memory..."

    # Use minimal threads + explicit memory cap for sort to prevent OOM
    # HT29_02N (150GB BAM) OOM'd twice: once at 32 threads, once at 8 threads
    # Problem: samtools sort default 768MB/thread + merge buffers exceed 60GB
    # Fix: 4 threads × 256MB = 1GB sort buffer — slower but safe
    SORT_THREADS=4
    SORT_MEM="256M"

    echo "  Sort config: $SORT_THREADS threads, $SORT_MEM per thread"

    samtools collate -@ $SORT_THREADS -o "$WORKDIR/bams/collate.bam" "$INPUT_BAM"
    rm -f "$INPUT_BAM" "${INPUT_BAM}.bai"
    echo "  collate done. Freed input BAM."

    samtools fixmate -@ $SORT_THREADS -m "$WORKDIR/bams/collate.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/collate.bam"
    echo "  fixmate done."

    samtools sort -m $SORT_MEM -@ $SORT_THREADS -o "$WORKDIR/bams/fixmate.sorted.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/fixmate.bam"
    echo "  sort done."

    samtools markdup -@ $SORT_THREADS -s "$WORKDIR/bams/fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt"
    rm -f "$WORKDIR/bams/fixmate.sorted.bam"
    echo "  markdup done."

    samtools index -@ $THREADS "$MARKDUP_BAM"
    rm -f "$INPUT_BAM" "${INPUT_BAM}.bai"

    # Upload markdup BAM (useful for P2/P3 reuse)
    echo "  Uploading markdup BAM..."
    aws s3 cp --quiet "$MARKDUP_BAM" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.markdup.bam"
    aws s3 cp --quiet "${MARKDUP_BAM}.bai" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.markdup.bam.bai"
    aws s3 cp --quiet "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.markdup_metrics.txt"
    echo "MarkDup complete."
else
    echo "[3] Skipping markdup (BAM_TYPE=$BAM_TYPE)"
    mv "$INPUT_BAM" "$MARKDUP_BAM"
    [ -f "${INPUT_BAM}.bai" ] && mv "${INPUT_BAM}.bai" "${MARKDUP_BAM}.bai" || samtools index -@ $THREADS "$MARKDUP_BAM"
fi

# ---- Calmd ----
CALMD_BAM="$WORKDIR/bams/${SAMPLE}.calmd.bam"
echo "[4] Running samtools calmd..."
samtools calmd -b -@ $THREADS "$MARKDUP_BAM" "$REF" > "$CALMD_BAM" 2> /dev/null
samtools index -@ $THREADS "$CALMD_BAM"
rm -f "$MARKDUP_BAM" "${MARKDUP_BAM}.bai"
echo "Calmd done."

# ---- Revelio ----
REVELIO_BAM="$WORKDIR/bams/${SAMPLE}.revelio.sorted.bam"
echo "[5] Running Revelio (12 chunks parallel)..."

TOTAL_READS=$(samtools view -c "$CALMD_BAM")
NUM_CHUNKS=12
CHUNK_SIZE=$(( (TOTAL_READS + NUM_CHUNKS - 1) / NUM_CHUNKS ))
THREADS_PER_CHUNK=3

echo "  Total reads: $TOTAL_READS, Chunks: $NUM_CHUNKS, Reads/chunk: $CHUNK_SIZE"

REVELIO_DIR="$WORKDIR/revelio_tmp"
mkdir -p "$REVELIO_DIR"
samtools view -H "$CALMD_BAM" > "$REVELIO_DIR/header.sam"

# Split, process, collect
samtools view "$CALMD_BAM" | split -l $CHUNK_SIZE -d -a 3 - "$REVELIO_DIR/chunk_"

REVELIO_OUTPUTS=()
for chunk in "$REVELIO_DIR"/chunk_*; do
    chunk_bam="${chunk}.bam"
    cat "$REVELIO_DIR/header.sam" "$chunk" | samtools view -b -o "$chunk_bam" -
    samtools index "$chunk_bam"
    rm -f "$chunk"

    output_bam="${chunk}.revelio.bam"
    REVELIO_OUTPUTS+=("$output_bam")
    python3 /opt/revelio.py -T $THREADS_PER_CHUNK -Q "$chunk_bam" "$output_bam" &
done

wait
echo "  All chunks done. Merging..."

samtools cat -o "$REVELIO_DIR/merged.bam" "${REVELIO_OUTPUTS[@]}"
samtools sort -@ $THREADS -o "$REVELIO_BAM" "$REVELIO_DIR/merged.bam"
samtools index -@ $THREADS "$REVELIO_BAM"
rm -rf "$REVELIO_DIR"
rm -f "$CALMD_BAM" "${CALMD_BAM}.bai"
echo "Revelio done."

# ---- BCFtools variant calling ----
VCF="$WORKDIR/results/${SAMPLE}.bcftools.vcf"
echo "[6] Running BCFtools mpileup + call..."
bcftools mpileup \
    --threads $THREADS \
    -Ou \
    -q 20 \
    -Q 20 \
    -f "$REF" \
    "$REVELIO_BAM" \
| bcftools call \
    --threads $THREADS \
    -mv \
    -Ov \
    -o "$VCF"

bgzip -c "$VCF" > "${VCF}.gz"
tabix -p vcf "${VCF}.gz"
VARIANT_COUNT=$(grep -vc '^#' "$VCF" || echo "0")
echo "Variant calling done: $VARIANT_COUNT variants"

# ---- Upload results ----
echo ""
echo "Uploading results to S3..."
aws s3 cp --quiet "$VCF" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.bcftools.vcf"
aws s3 cp --quiet "${VCF}.gz" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.bcftools.vcf.gz"
aws s3 cp --quiet "${VCF}.gz.tbi" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.bcftools.vcf.gz.tbi"
aws s3 cp --quiet "$REVELIO_BAM" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.revelio.sorted.bam"
aws s3 cp --quiet "${REVELIO_BAM}.bai" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.revelio.sorted.bam.bai"

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE ($VARIANT_COUNT variants)"
echo "=========================================="
echo "Finished: $(date)"
echo "Results: $S3_OUTPUT/p1_dna_snp/${SAMPLE}/"
