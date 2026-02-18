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
    echo "[3] Running samtools markdup (5-way chromosome split, 2-stage)..."

    # OOM history (HT29_02N 140GB BAM, 2B reads, 60GB container):
    #   v3-v5: Various pipe approaches OOM — hash table too large for full/half genome
    #   v6: 2-way split pipe OOM — chr1-11 = 63%, hash ~45GB + page cache
    #   v7: 3-way split pipe OOM — Group A (chr1-6, 60GB) still too big; sort temp files
    #        + markdup hash + page cache coexist in pipe → 15 min after sort merge = OOM
    # Fix (v8): 5-way split + TWO-STAGE processing.
    #   Stage 1: collate|fixmate|sort → file (no markdup in memory)
    #   Stage 2: markdup on coordinate-sorted file (streaming, low memory)
    #   This ensures sort temp files and markdup hash never coexist in memory.
    #   Groups: A=chr1-2(~16%), B=chr3-5(~19%), C=chr6-8-9(~20%), D=chr10-14(~20%), E=chr15-Y(~25%)
    #   Max group BAM ~35GB, max hash ~18GB — each stage independently fits in 60GB.
    PIPE_THREADS=2
    SORT_MEM="256M"
    GROUP_A="chr1 chr2"
    GROUP_B="chr3 chr4 chr5"
    GROUP_C="chr6 chr7 chr8 chr9"
    GROUP_D="chr10 chr11 chr12 chr13 chr14"
    GROUP_E="chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

    mkdir -p "$WORKDIR/bams/tmp"

    # Index needed for region extraction
    if [ ! -f "${INPUT_BAM}.bai" ]; then
        echo "  Indexing input BAM..."
        samtools index -@ $THREADS "$INPUT_BAM"
    fi

    # Step 1: Extract each group to a separate BAM file
    for GNAME in A B C D E; do
        eval "REGIONS=\$GROUP_${GNAME}"
        echo "  Extracting group ${GNAME} (${REGIONS})..."
        samtools view -b -@ 4 "$INPUT_BAM" $REGIONS > "$WORKDIR/bams/group${GNAME}.bam"
    done

    echo "  Group sizes:"
    ls -lh "$WORKDIR/bams/group"*.bam

    # Step 2: DELETE original BAM to free 140GB from page cache
    echo "  Deleting original BAM to free memory..."
    rm -f "$INPUT_BAM" "${INPUT_BAM}.bai"
    sync
    echo 3 > /proc/sys/vm/drop_caches 2>/dev/null || true

    # Step 3: Process each group with 2-stage approach
    for GNAME in A B C D E; do
        GROUP_BAM="$WORKDIR/bams/group${GNAME}.bam"
        SORTED_BAM="$WORKDIR/bams/group${GNAME}.sorted.bam"
        OUTPUT_BAM="$WORKDIR/bams/group${GNAME}.markdup.bam"
        GROUP_SIZE=$(du -h "$GROUP_BAM" | cut -f1)
        echo "  Group ${GNAME}: processing (${GROUP_SIZE})..."

        # Stage 1: collate + fixmate + sort → coordinate-sorted file
        # No markdup in memory during this stage
        echo "    Stage 1: collate|fixmate|sort..."
        samtools collate -u -@ $PIPE_THREADS -O "$GROUP_BAM" "$WORKDIR/bams/tmp/collate_${GNAME}" | \
            samtools fixmate -u -@ $PIPE_THREADS -m - - | \
            samtools sort -m $SORT_MEM -@ $PIPE_THREADS -T "$WORKDIR/bams/tmp/sort_${GNAME}" \
            -o "$SORTED_BAM" -

        # Free memory: delete input BAM + drop page cache before markdup
        rm -f "$GROUP_BAM"
        rm -f "$WORKDIR/bams/tmp/"*
        sync
        echo 3 > /proc/sys/vm/drop_caches 2>/dev/null || true

        # Stage 2: markdup on coordinate-sorted input (streaming, low memory)
        echo "    Stage 2: markdup..."
        samtools markdup -@ $PIPE_THREADS -s "$SORTED_BAM" "$OUTPUT_BAM" \
            2>> "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt"

        rm -f "$SORTED_BAM"
        sync
        echo 3 > /proc/sys/vm/drop_caches 2>/dev/null || true
        echo "  Group ${GNAME} done."
    done

    # Step 4: Merge groups
    echo "  Merging 5 groups..."
    samtools merge -@ $THREADS "$MARKDUP_BAM" \
        "$WORKDIR/bams/groupA.markdup.bam" \
        "$WORKDIR/bams/groupB.markdup.bam" \
        "$WORKDIR/bams/groupC.markdup.bam" \
        "$WORKDIR/bams/groupD.markdup.bam" \
        "$WORKDIR/bams/groupE.markdup.bam"
    rm -f "$WORKDIR/bams/group"*.markdup.bam
    rm -rf "$WORKDIR/bams/tmp"

    samtools index -@ $THREADS "$MARKDUP_BAM"
    echo "  MarkDup done."

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
