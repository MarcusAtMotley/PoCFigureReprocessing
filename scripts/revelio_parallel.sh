#!/bin/bash
# Parallel Revelio - splits BAM by read count, processes chunks in parallel, merges
#
# Usage: revelio_parallel.sh -i INPUT_BAM -o OUTPUT_BAM [-n NUM_CHUNKS] [-t THREADS_PER_CHUNK]
#
# This bypasses Python's GIL by running multiple revelio processes

set -e

# Defaults
NUM_CHUNKS=12
THREADS_PER_CHUNK=3
REVELIO_SCRIPT=$(dirname "$0")/../pipelines/modules/local/revelio/bin/revelio.py

usage() {
    echo "Usage: $0 -i INPUT_BAM -o OUTPUT_BAM [-n NUM_CHUNKS] [-t THREADS_PER_CHUNK]"
    echo ""
    echo "Options:"
    echo "  -i  Input BAM file (must have MD tags from samtools calmd)"
    echo "  -o  Output BAM file"
    echo "  -n  Number of chunks to split into (default: 12)"
    echo "  -t  Threads per revelio chunk (default: 3)"
    echo "  -r  Path to revelio.py script (auto-detected by default)"
    exit 1
}

while getopts "i:o:n:t:r:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        o) OUTPUT_BAM="$OPTARG" ;;
        n) NUM_CHUNKS="$OPTARG" ;;
        t) THREADS_PER_CHUNK="$OPTARG" ;;
        r) REVELIO_SCRIPT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "$INPUT_BAM" ] || [ -z "$OUTPUT_BAM" ]; then
    echo "ERROR: Input and output BAM are required"
    usage
fi

if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: Input BAM not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "$REVELIO_SCRIPT" ]; then
    echo "ERROR: Revelio script not found: $REVELIO_SCRIPT"
    exit 1
fi

# Create temp directory for chunks
WORK_DIR=$(mktemp -d -p /data/p5_work revelio_parallel.XXXXXX)
trap "rm -rf $WORK_DIR" EXIT

echo "=========================================="
echo "Parallel Revelio"
echo "=========================================="
echo "Input:  $INPUT_BAM"
echo "Output: $OUTPUT_BAM"
echo "Chunks: $NUM_CHUNKS"
echo "Threads per chunk: $THREADS_PER_CHUNK"
echo "Work dir: $WORK_DIR"
echo ""

# Step 1: Count reads and calculate chunk size
echo "[1/4] Counting reads..."
TOTAL_READS=$(samtools view -c "$INPUT_BAM")
CHUNK_SIZE=$(( (TOTAL_READS + NUM_CHUNKS - 1) / NUM_CHUNKS ))
echo "Total reads: $TOTAL_READS"
echo "Chunk size: ~$CHUNK_SIZE reads"
echo ""

# Step 2: Extract header
echo "[2/4] Splitting BAM into $NUM_CHUNKS chunks..."
samtools view -H "$INPUT_BAM" > "$WORK_DIR/header.sam"

# Split BAM into chunks by streaming through split
# Using SAM format for splitting, then convert back to BAM
samtools view "$INPUT_BAM" | split -l "$CHUNK_SIZE" -d -a 3 - "$WORK_DIR/chunk_"

# Convert each chunk to BAM with header and index
CHUNK_FILES=()
for chunk in "$WORK_DIR"/chunk_*; do
    chunk_name=$(basename "$chunk")
    chunk_bam="$WORK_DIR/${chunk_name}.bam"
    cat "$WORK_DIR/header.sam" "$chunk" | samtools view -b -o "$chunk_bam" -
    samtools index "$chunk_bam"  # Index for revelio
    rm "$chunk"  # Clean up SAM chunk
    CHUNK_FILES+=("$chunk_bam")
done
echo "Created and indexed ${#CHUNK_FILES[@]} chunks"
echo ""

# Step 3: Run revelio on each chunk in parallel
echo "[3/4] Running Revelio on ${#CHUNK_FILES[@]} chunks in parallel..."
REVELIO_PIDS=()
for chunk_bam in "${CHUNK_FILES[@]}"; do
    chunk_name=$(basename "$chunk_bam" .bam)
    output_bam="$WORK_DIR/${chunk_name}.revelio.bam"

    python3 "$REVELIO_SCRIPT" -T "$THREADS_PER_CHUNK" -Q "$chunk_bam" "$output_bam" &
    REVELIO_PIDS+=($!)
    echo "  Started revelio on $chunk_name (PID: $!)"
done

# Wait for all revelio processes
echo ""
echo "Waiting for ${#REVELIO_PIDS[@]} revelio processes..."
FAILED=0
for pid in "${REVELIO_PIDS[@]}"; do
    if ! wait "$pid"; then
        echo "ERROR: Revelio process $pid failed"
        FAILED=1
    fi
done

if [ "$FAILED" -eq 1 ]; then
    echo "ERROR: One or more revelio processes failed"
    exit 1
fi
echo "All revelio processes completed"
echo ""

# Step 4: Merge results and sort
echo "[4/4] Merging and sorting results..."
REVELIO_OUTPUTS=()
for chunk_bam in "${CHUNK_FILES[@]}"; do
    chunk_name=$(basename "$chunk_bam" .bam)
    REVELIO_OUTPUTS+=("$WORK_DIR/${chunk_name}.revelio.bam")
done

# Concatenate (not merge - preserves order) and sort
samtools cat -o "$WORK_DIR/merged.bam" "${REVELIO_OUTPUTS[@]}"
samtools sort -@ 8 -o "$OUTPUT_BAM" "$WORK_DIR/merged.bam"
samtools index -@ 8 "$OUTPUT_BAM"

echo ""
echo "=========================================="
echo "Parallel Revelio complete"
echo "=========================================="
echo "Output: $OUTPUT_BAM"
echo "Index:  ${OUTPUT_BAM}.bai"
