#!/bin/bash
# Self-contained DNA pipeline for AWS Batch
# Runs one sample end-to-end: download → align → markdup → calmd → revelio → bcftools → upload
#
# Usage: batch_dna_pipeline.sh <SAMPLE> <SINGLE_END> <FASTQ_1_PATHS> [FASTQ_2_PATHS]
#   SAMPLE: sample name
#   SINGLE_END: "true" | "false"
#   FASTQ_1_PATHS: semicolon-separated S3 paths for R1 (or SE)
#   FASTQ_2_PATHS: semicolon-separated S3 paths for R2 (omit for SE)
#
# Example (paired-end, multi-lane):
#   batch_dna_pipeline.sh CoB_08R false "s3://bucket/R1_L001.fq;s3://bucket/R1_L002.fq" "s3://bucket/R2_L001.fq;s3://bucket/R2_L002.fq"
#
# Example (single-end):
#   batch_dna_pipeline.sh HT29_21W true "s3://bucket/SE.fq"
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE SINGLE_END FASTQ_1 [FASTQ_2]}"
SINGLE_END="${2:?Missing SINGLE_END (true/false)}"
FQ1_PATHS="${3:?Missing FASTQ_1 paths}"
FQ2_PATHS="${4:-}"

THREADS=$(nproc)
WORKDIR=/scratch
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
S3_REFS="s3://motleybio/Resources"
REF="$WORKDIR/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BISCUIT_INDEX="$WORKDIR/refs/biscuit_index/GRCh38_full_analysis_set_plus_decoy_hla.fa"

echo "=========================================="
echo "AWS Batch DNA Pipeline"
echo "=========================================="
echo "Sample: $SAMPLE"
echo "Single-end: $SINGLE_END"
echo "Threads: $THREADS"
echo "R1 paths: $FQ1_PATHS"
echo "R2 paths: $FQ2_PATHS"
echo "Started: $(date)"
echo ""

# ---- Step 1: Download references ----
echo "[1/7] Downloading references..."
mkdir -p "$WORKDIR/refs/biscuit_index" "$WORKDIR/fastq" "$WORKDIR/bams" "$WORKDIR/results"

if [ ! -f "$REF" ]; then
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" "$REF"
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" "${REF}.fai"
fi

if [ ! -f "${BISCUIT_INDEX}.bis.pac" ]; then
    aws s3 cp --quiet "$S3_REFS/biscuit_reference_genome/" "$WORKDIR/refs/biscuit_index/" --recursive
fi
echo "References ready."

# ---- Step 2: Download and merge FASTQs ----
echo "[2/7] Downloading FASTQs..."
FQ1="$WORKDIR/fastq/${SAMPLE}_R1.fastq"
FQ2="$WORKDIR/fastq/${SAMPLE}_R2.fastq"

# Download R1 files
IDX=0
IFS=';' read -ra R1_URLS <<< "$FQ1_PATHS"
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
if [ "$SINGLE_END" = "false" ] && [ -n "$FQ2_PATHS" ]; then
    IDX=0
    IFS=';' read -ra R2_URLS <<< "$FQ2_PATHS"
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

echo "FASTQs ready."
ls -lh "$WORKDIR/fastq/"

# ---- Step 3: Biscuit alignment ----
SORTED_BAM="$WORKDIR/bams/${SAMPLE}.sorted.bam"
echo "[3/7] Running Biscuit alignment..."
if [ "$SINGLE_END" = "true" ]; then
    biscuit align -@ $THREADS "$BISCUIT_INDEX" "$FQ1" \
        | samtools sort -@ 8 -o "$SORTED_BAM" -
else
    biscuit align -@ $THREADS "$BISCUIT_INDEX" "$FQ1" "$FQ2" \
        | samtools sort -@ 8 -o "$SORTED_BAM" -
fi
samtools index -@ $THREADS "$SORTED_BAM"
echo "Alignment done: $(ls -lh $SORTED_BAM)"

# Upload sorted BAM as checkpoint
echo "  Uploading alignment checkpoint..."
aws s3 cp --quiet "$SORTED_BAM" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.sorted.bam"
aws s3 cp --quiet "${SORTED_BAM}.bai" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.sorted.bam.bai"

# Free FASTQs
rm -f "$WORKDIR/fastq/"*

# ---- Step 4: Mark duplicates ----
MARKDUP_BAM="$WORKDIR/bams/${SAMPLE}.markdup.bam"
echo "[4/7] Running samtools markdup..."
samtools collate -@ $THREADS -o "$WORKDIR/bams/collate.bam" "$SORTED_BAM"
samtools fixmate -@ $THREADS -m "$WORKDIR/bams/collate.bam" "$WORKDIR/bams/fixmate.bam"
rm -f "$WORKDIR/bams/collate.bam"
samtools sort -@ $THREADS -o "$WORKDIR/bams/fixmate.sorted.bam" "$WORKDIR/bams/fixmate.bam"
rm -f "$WORKDIR/bams/fixmate.bam"
samtools markdup -@ $THREADS -s "$WORKDIR/bams/fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt"
rm -f "$WORKDIR/bams/fixmate.sorted.bam"
samtools index -@ $THREADS "$MARKDUP_BAM"
rm -f "$SORTED_BAM" "${SORTED_BAM}.bai"
echo "MarkDup done."

# ---- Step 5: Calmd (add MD tags) ----
CALMD_BAM="$WORKDIR/bams/${SAMPLE}.calmd.bam"
echo "[5/7] Running samtools calmd..."
samtools calmd -b -@ $THREADS "$MARKDUP_BAM" "$REF" > "$CALMD_BAM" 2> /dev/null
samtools index -@ $THREADS "$CALMD_BAM"
rm -f "$MARKDUP_BAM" "${MARKDUP_BAM}.bai"
echo "Calmd done."

# ---- Step 6: Revelio ----
REVELIO_BAM="$WORKDIR/bams/${SAMPLE}.revelio.sorted.bam"
echo "[6/7] Running Revelio (12 chunks parallel)..."

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

# ---- Step 7: BCFtools variant calling ----
VCF="$WORKDIR/results/${SAMPLE}.bcftools.vcf"
echo "[7/7] Running BCFtools mpileup + call..."
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
aws s3 cp --quiet "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt" "$S3_OUTPUT/p1_dna_snp/${SAMPLE}/${SAMPLE}.markdup_metrics.txt"

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE ($VARIANT_COUNT variants)"
echo "=========================================="
echo "Finished: $(date)"
echo "Results: $S3_OUTPUT/p1_dna_snp/${SAMPLE}/"
