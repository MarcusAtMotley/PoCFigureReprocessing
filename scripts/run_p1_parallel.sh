#!/bin/bash
# Run P1 DNA SNP Pipeline - 3 samples in parallel at 16 threads each
set -e

SCRIPT_DIR=$(dirname "$0")
THREADS_PER_SAMPLE=16
WORKDIR=/data/dna_work
RESULTS=/data/results
FASTQ_DIR=/data/fastq
S3_BASE="s3://motleybio/Datasets/GOOGLE_CLOUD/run_motley26/results/rna_deconvolution/cutadapt"

# Source the pipeline functions
source "$SCRIPT_DIR/simple_dna_pipeline.sh"

# Override thread count
THREADS=$THREADS_PER_SAMPLE

mkdir -p "$FASTQ_DIR" "$WORKDIR" "$RESULTS/p1_dna_snp"

echo "=========================================="
echo "P1 DNA SNP - Parallel (3 samples @ ${THREADS_PER_SAMPLE} threads)"
echo "=========================================="

# Function to download and merge multi-lane FASTQs
download_and_merge() {
    local SAMPLE=$1
    local S3_PATTERN=$2
    local OUT_R1="$FASTQ_DIR/${SAMPLE}_merged_R1.fastq"
    local OUT_R2="$FASTQ_DIR/${SAMPLE}_merged_R2.fastq"

    if [ -f "$OUT_R1" ] && [ -f "$OUT_R2" ]; then
        echo "[$SAMPLE] Merged FASTQs exist, skipping download..."
        return 0
    fi

    echo "[$SAMPLE] Downloading from S3..."
    local TMPDIR="$FASTQ_DIR/tmp_${SAMPLE}"
    mkdir -p "$TMPDIR"

    # Download all lanes
    aws s3 cp "$S3_BASE/" "$TMPDIR/" --recursive --exclude "*" --include "${S3_PATTERN}*unbarcoded*.fastq" --quiet

    echo "[$SAMPLE] Merging lanes..."
    cat "$TMPDIR"/*_R1.cutadapt.fastq > "$OUT_R1"
    cat "$TMPDIR"/*_R2.cutadapt.fastq > "$OUT_R2"

    rm -rf "$TMPDIR"
    echo "[$SAMPLE] Download complete: $(ls -lh $OUT_R1 | awk '{print $5}')"
}

# Function to process a single sample
process_sample_wrapper() {
    local SAMPLE=$1
    local FQ1=$2
    local FQ2=$3
    local LOGFILE="$WORKDIR/${SAMPLE}.log"

    echo "[$SAMPLE] Starting processing..."
    process_dna_sample "$SAMPLE" "$FQ1" "$FQ2" "false" > "$LOGFILE" 2>&1
    echo "[$SAMPLE] Complete!"
}

# Sample 1: CoB_08L_3A2_DNA-EM (already local)
SAMPLE1="CoB_08L_3A2_DNA-EM"
FQ1_1="$FASTQ_DIR/CoB_08L_3A2_DNA-EM_S1_merged_R1.fastq"
FQ1_2="$FASTQ_DIR/CoB_08L_3A2_DNA-EM_S1_merged_R2.fastq"

# Sample 2: CoM_08M_3A2_DNA-EM (download)
SAMPLE2="CoM_08M_3A2_DNA-EM"
FQ2_1="$FASTQ_DIR/${SAMPLE2}_merged_R1.fastq"
FQ2_2="$FASTQ_DIR/${SAMPLE2}_merged_R2.fastq"

# Sample 3: CoB_08R_3A2_TNA-mRT-EM (download) - DNA channel
SAMPLE3="CoB_08R_3A2_TNA-mRT-EM"
FQ3_1="$FASTQ_DIR/${SAMPLE3}_merged_R1.fastq"
FQ3_2="$FASTQ_DIR/${SAMPLE3}_merged_R2.fastq"

echo ""
echo "Samples to process:"
echo "  1. $SAMPLE1 (local)"
echo "  2. $SAMPLE2 (downloading)"
echo "  3. $SAMPLE3 (downloading)"
echo ""

# Start downloads in background
download_and_merge "$SAMPLE2" "CoM_08M_3A2_DNA-EM_S2" &
DL_PID2=$!

download_and_merge "$SAMPLE3" "CoB_08R_3A2_TNA-mRT-EM_S3" &
DL_PID3=$!

# Start sample 1 immediately (already local)
echo ""
echo "=== Starting $SAMPLE1 (local) ==="
process_sample_wrapper "$SAMPLE1" "$FQ1_1" "$FQ1_2" &
PROC_PID1=$!

# Wait for sample 2 download, then start processing
wait $DL_PID2
echo ""
echo "=== Starting $SAMPLE2 ==="
process_sample_wrapper "$SAMPLE2" "$FQ2_1" "$FQ2_2" &
PROC_PID2=$!

# Wait for sample 3 download, then start processing
wait $DL_PID3
echo ""
echo "=== Starting $SAMPLE3 ==="
process_sample_wrapper "$SAMPLE3" "$FQ3_1" "$FQ3_2" &
PROC_PID3=$!

# Wait for all processing to complete
echo ""
echo "All samples running in parallel. Waiting for completion..."
wait $PROC_PID1 && echo "✓ $SAMPLE1 done"
wait $PROC_PID2 && echo "✓ $SAMPLE2 done"
wait $PROC_PID3 && echo "✓ $SAMPLE3 done"

echo ""
echo "=========================================="
echo "P1 DNA SNP Complete!"
echo "=========================================="
echo "Results:"
ls -la "$RESULTS/p1_dna_snp/"
