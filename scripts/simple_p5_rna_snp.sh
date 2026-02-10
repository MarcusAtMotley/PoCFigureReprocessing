#!/bin/bash
# Simple P5 RNA SNP Pipeline - processes P4 BAMs through Revelio + BCFtools
# No Nextflow, just straightforward bash
#
# Key changes from original:
# - BCFtools replaces LoFreq (50-100x faster)
# - Conditional Revelio (only for bisulfite/EM treated samples)
# - Downsampling support for SingleAnalyte samples

set -e

# Configuration
THREADS=40  # Use most CPUs, leave some for system/I/O
WORKDIR=/data/p5_work
RESULTS=/data/results/p5_rna_snp
REF=~/references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa
S3_INPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p4_rna_counts"
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p5_rna_snp"
REVELIO_SCRIPT=$(dirname "$0")/../pipelines/modules/local/revelio/bin/revelio.py
REVELIO_PARALLEL=$(dirname "$0")/revelio_parallel.sh
DOWNSAMPLE_SCRIPT=$(dirname "$0")/downsample_bam.sh

# Downsampling targets
DOWNSAMPLE_RNA=30000000  # 30M reads for SingleAnalyte RNA to match mTNA

# Create directories
mkdir -p "$WORKDIR" "$RESULTS"

echo "=========================================="
echo "P5 RNA SNP Pipeline (BCFtools)"
echo "=========================================="
echo "Threads: $THREADS"
echo "Work dir: $WORKDIR"
echo "Results: $RESULTS"
echo ""

# Check reference exists
if [ ! -f "$REF" ]; then
    echo "ERROR: Reference not found at $REF"
    exit 1
fi

# Index reference if needed
if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference..."
    samtools faidx "$REF"
fi

# All samples get Revelio for pipeline consistency
# Even non-bisulfite samples - ensures identical processing so differences reflect biology
needs_revelio() {
    # Always return true - apply Revelio to ALL samples for comparability
    return 0
}

# Determine if sample needs downsampling
needs_downsampling() {
    local SAMPLE=$1
    case "$SAMPLE" in
        *"1A3"*)
            # SingleAnalyte samples (1A3 in name) need downsampling
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

process_sample() {
    local SAMPLE=$1
    echo ""
    echo "=========================================="
    echo "Processing: $SAMPLE"
    echo "=========================================="

    local SAMPLE_DIR="$WORKDIR/$SAMPLE"
    local OUT_DIR="$RESULTS/$SAMPLE"
    mkdir -p "$SAMPLE_DIR" "$OUT_DIR"

    local INPUT_BAM="$SAMPLE_DIR/${SAMPLE}.bam"
    local WORKING_BAM="$INPUT_BAM"  # May change if downsampled
    local CALMD_BAM="$SAMPLE_DIR/${SAMPLE}.calmd.bam"
    local REVELIO_BAM="$SAMPLE_DIR/${SAMPLE}.revelio.bam"
    local SORTED_BAM="$SAMPLE_DIR/${SAMPLE}.revelio.sorted.bam"
    local VCF="$OUT_DIR/${SAMPLE}.bcftools.vcf"

    # Step 1: Download BAM from S3 (if not present)
    if [ ! -f "$INPUT_BAM" ]; then
        echo "[1/6] Downloading BAM from S3..."
        aws s3 cp --quiet "${S3_INPUT}/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam" "$INPUT_BAM"
    else
        echo "[1/6] BAM already downloaded, skipping..."
    fi

    # Step 2: Index input BAM
    if [ ! -f "${INPUT_BAM}.bai" ]; then
        echo "[2/6] Indexing input BAM..."
        samtools index -@ $THREADS "$INPUT_BAM"
    else
        echo "[2/6] BAM index exists, skipping..."
    fi

    # Step 3: Downsample if needed (SingleAnalyte samples)
    if needs_downsampling "$SAMPLE"; then
        local DOWNSAMPLED_BAM="$SAMPLE_DIR/${SAMPLE}.downsampled.bam"
        if [ ! -f "$DOWNSAMPLED_BAM" ]; then
            echo "[3/6] Downsampling to ${DOWNSAMPLE_RNA} reads..."
            "$DOWNSAMPLE_SCRIPT" -i "$INPUT_BAM" -o "$DOWNSAMPLED_BAM" -t "$DOWNSAMPLE_RNA" -@ "$THREADS"
        else
            echo "[3/6] Downsampled BAM exists, skipping..."
        fi
        WORKING_BAM="$DOWNSAMPLED_BAM"
    else
        echo "[3/6] No downsampling needed for TrinitySeq sample"
    fi

    # Step 4: Add MD tags with samtools calmd
    if [ ! -f "$CALMD_BAM" ]; then
        echo "[4/6] Running samtools calmd..."
        samtools calmd -b -@ $THREADS "$WORKING_BAM" "$REF" > "$CALMD_BAM" 2> "$SAMPLE_DIR/calmd.log"
        samtools index -@ $THREADS "$CALMD_BAM"
    else
        echo "[4/6] calmd BAM exists, skipping..."
    fi

    # Step 5: Run Revelio on ALL samples for pipeline consistency (PARALLEL)
    local FINAL_BAM="$SORTED_BAM"
    if [ ! -f "$SORTED_BAM" ]; then
        echo "[5/6] Running Parallel Revelio (12 chunks, applied to ALL samples)..."
        "$REVELIO_PARALLEL" -i "$CALMD_BAM" -o "$SORTED_BAM" -n 12 -t 3
    else
        echo "[5/6] Revelio BAM exists, skipping..."
    fi

    # Step 6: Run BCFtools for variant calling
    # RNA SNP uses lower quality thresholds: -q 13 -Q 13
    if [ ! -f "$VCF" ]; then
        echo "[6/6] Running BCFtools mpileup + call..."
        bcftools mpileup \
            --threads $THREADS \
            -Ou \
            -q 13 \
            -Q 13 \
            -f "$REF" \
            "$FINAL_BAM" \
        | bcftools call \
            --threads $THREADS \
            -mv \
            -Ov \
            -o "$VCF"

        # Compress and index VCF
        bgzip -c "$VCF" > "${VCF}.gz"
        tabix -p vcf "${VCF}.gz"
    else
        echo "[6/6] VCF exists, skipping..."
    fi

    echo "Done: $SAMPLE -> $VCF"

    # Optional: Upload results to S3
    # aws s3 cp "$VCF" "${S3_OUTPUT}/${SAMPLE}/"
    # aws s3 cp "${VCF}.gz" "${S3_OUTPUT}/${SAMPLE}/"
}

# TrinitySeq RNA Samples (6 total - only samples with P4 BAMs on S3)
# CoB: 1 mTNA
# CoM: 1 mTNA
# HT29: 4 HairyTNA
TRINITYSEQ_SAMPLES=(
    # CoB mTNA
    "CoB_08R_3A2_TNA-mRT-EM"
    # CoM mTNA
    "CoM_08S_3A2_TNA-mRT-EM"
    # HT29 HairyTNA
    "HT29_21S_3A2_bsTNA-HP"
    "HT29_21T_3A2_bsTNA-HP"
    "HT29_21U_3A2_bsRNA-HP"
    "HT29_21V_3A2_bsRNA-HP"
)

# SingleAnalyte RNA Samples (3 total) - will be downsampled
SINGLEANALYTE_SAMPLES=(
    "CoB_02W_1A3_1RNA"
    "CoM_02V_1A3_1RNA"
    "HT29_02T_1A3_1RNA"
)

# Default: process all samples
ALL_SAMPLES=("${TRINITYSEQ_SAMPLES[@]}" "${SINGLEANALYTE_SAMPLES[@]}")

# Parse command line arguments
SAMPLES_TO_PROCESS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --trinityseq-only)
            SAMPLES_TO_PROCESS=("${TRINITYSEQ_SAMPLES[@]}")
            shift
            ;;
        --singleanalyte-only)
            SAMPLES_TO_PROCESS=("${SINGLEANALYTE_SAMPLES[@]}")
            shift
            ;;
        --sample)
            SAMPLES_TO_PROCESS+=("$2")
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --trinityseq-only    Process only TrinitySeq samples (12)"
            echo "  --singleanalyte-only Process only SingleAnalyte samples (3)"
            echo "  --sample SAMPLE      Process specific sample (can be repeated)"
            echo "  --help               Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# If no options specified, process all
if [ ${#SAMPLES_TO_PROCESS[@]} -eq 0 ]; then
    SAMPLES_TO_PROCESS=("${ALL_SAMPLES[@]}")
fi

echo "Samples to process: ${#SAMPLES_TO_PROCESS[@]}"
echo ""

# Process each sample
for SAMPLE in "${SAMPLES_TO_PROCESS[@]}"; do
    process_sample "$SAMPLE"
done

echo ""
echo "=========================================="
echo "P5 Complete!"
echo "=========================================="
echo "Results in: $RESULTS"
ls -la "$RESULTS"
