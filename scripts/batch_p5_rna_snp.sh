#!/bin/bash
# P5 RNA SNP pipeline for AWS Batch — from STAR-aligned BAM
# Downloads P4 BAM, runs calmd → revelio → bcftools mpileup+call
#
# Usage: batch_p5_rna_snp.sh <SAMPLE> <BAM_S3_PATH>
#   SAMPLE:      sample name
#   BAM_S3_PATH: S3 path to STAR-aligned BAM (from P4)
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE BAM_S3_PATH}"
BAM_S3="${2:?Missing BAM_S3_PATH}"

THREADS=$(nproc)
WORKDIR=/scratch
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
S3_REFS="s3://motleybio/Resources"
REF="$WORKDIR/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa"

echo "=========================================="
echo "P5 RNA SNP — from STAR BAM"
echo "=========================================="
echo "Sample: $SAMPLE"
echo "BAM source: $BAM_S3"
echo "Threads: $THREADS"
echo "Started: $(date)"
echo ""

# ---- Clean scratch ----
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
echo "[2] Downloading STAR BAM..."
INPUT_BAM="$WORKDIR/bams/${SAMPLE}.bam"
aws s3 cp --quiet "$BAM_S3" "$INPUT_BAM"
samtools index -@ $THREADS "$INPUT_BAM"
echo "BAM ready: $(ls -lh $INPUT_BAM | awk '{print $5}')"

# ---- Calmd ----
CALMD_BAM="$WORKDIR/bams/${SAMPLE}.calmd.bam"
echo "[3] Running samtools calmd..."
samtools calmd -b -@ $THREADS "$INPUT_BAM" "$REF" > "$CALMD_BAM" 2>/dev/null
samtools index -@ $THREADS "$CALMD_BAM"
rm -f "$INPUT_BAM" "${INPUT_BAM}.bai"
echo "Calmd done."

# ---- Revelio (chromosome-based splitting) ----
REVELIO_BAM="$WORKDIR/bams/${SAMPLE}.revelio.sorted.bam"
echo "[4] Running Revelio (chromosome groups)..."

REVELIO_DIR="$WORKDIR/revelio_tmp"
mkdir -p "$REVELIO_DIR"

# 6 chromosome groups for parallel revelio
GROUP_A_REGIONS="chr1 chr2 chr3"
GROUP_B_REGIONS="chr4 chr5 chr6 chr7"
GROUP_C_REGIONS="chr8 chr9 chr10 chr11"
GROUP_D_REGIONS="chr12 chr13 chr14 chr15 chr16"
GROUP_E_REGIONS="chr17 chr18 chr19 chr20 chr21 chr22"
GROUP_F_REGIONS="chrX chrY chrM"

THREADS_PER_GROUP=$(( THREADS / 6 ))
[ "$THREADS_PER_GROUP" -lt 1 ] && THREADS_PER_GROUP=1

for group in A B C D E F; do
    REGIONS_VAR="GROUP_${group}_REGIONS"
    REGIONS="${!REGIONS_VAR}"
    GROUP_BAM="$REVELIO_DIR/group_${group}.bam"
    OUTPUT_BAM="$REVELIO_DIR/group_${group}.revelio.bam"

    samtools view -b -@ 4 -o "$GROUP_BAM" "$CALMD_BAM" $REGIONS
    samtools index -@ 4 "$GROUP_BAM"
    echo "  Group $group: started revelio"
    python3 /opt/revelio.py -T $THREADS_PER_GROUP -Q "$GROUP_BAM" "$OUTPUT_BAM" &
done

wait
echo "  All groups done. Merging..."

# Merge revelio outputs
samtools cat -o "$REVELIO_DIR/merged.bam" \
    "$REVELIO_DIR/group_A.revelio.bam" \
    "$REVELIO_DIR/group_B.revelio.bam" \
    "$REVELIO_DIR/group_C.revelio.bam" \
    "$REVELIO_DIR/group_D.revelio.bam" \
    "$REVELIO_DIR/group_E.revelio.bam" \
    "$REVELIO_DIR/group_F.revelio.bam"

samtools sort -@ $THREADS -o "$REVELIO_BAM" "$REVELIO_DIR/merged.bam"
samtools index -@ $THREADS "$REVELIO_BAM"
rm -rf "$REVELIO_DIR"
rm -f "$CALMD_BAM" "${CALMD_BAM}.bai"
echo "Revelio done."

# ---- BCFtools variant calling ----
VCF="$WORKDIR/results/${SAMPLE}.bcftools.vcf"
echo "[5] Running BCFtools mpileup + call..."
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
aws s3 cp --quiet "$VCF" "$S3_OUTPUT/p5_rna_snp/${SAMPLE}/${SAMPLE}.bcftools.vcf"
aws s3 cp --quiet "${VCF}.gz" "$S3_OUTPUT/p5_rna_snp/${SAMPLE}/${SAMPLE}.bcftools.vcf.gz"
aws s3 cp --quiet "${VCF}.gz.tbi" "$S3_OUTPUT/p5_rna_snp/${SAMPLE}/${SAMPLE}.bcftools.vcf.gz.tbi"
aws s3 cp --quiet "$REVELIO_BAM" "$S3_OUTPUT/p5_rna_snp/${SAMPLE}/${SAMPLE}.revelio.sorted.bam"
aws s3 cp --quiet "${REVELIO_BAM}.bai" "$S3_OUTPUT/p5_rna_snp/${SAMPLE}/${SAMPLE}.revelio.sorted.bam.bai"

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE ($VARIANT_COUNT variants)"
echo "=========================================="
echo "Finished: $(date)"
echo "Results: $S3_OUTPUT/p5_rna_snp/${SAMPLE}/"
