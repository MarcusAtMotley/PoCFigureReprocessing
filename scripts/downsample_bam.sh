#!/bin/bash
# Downsample BAM file to target number of reads
# Uses samtools view -s with seed-based random sampling

set -e

usage() {
    echo "Usage: $0 -i INPUT_BAM -o OUTPUT_BAM -t TARGET_READS [-s SEED] [-@ THREADS]"
    echo ""
    echo "Options:"
    echo "  -i  Input BAM file"
    echo "  -o  Output BAM file"
    echo "  -t  Target number of reads"
    echo "  -s  Random seed (default: 42)"
    echo "  -@  Number of threads (default: 8)"
    echo ""
    echo "Example:"
    echo "  $0 -i sample.bam -o sample.downsampled.bam -t 30000000"
    exit 1
}

# Defaults
SEED=42
THREADS=8

while getopts "i:o:t:s:@:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        o) OUTPUT_BAM="$OPTARG" ;;
        t) TARGET_READS="$OPTARG" ;;
        s) SEED="$OPTARG" ;;
        @) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_BAM" ] || [ -z "$OUTPUT_BAM" ] || [ -z "$TARGET_READS" ]; then
    echo "ERROR: Missing required arguments"
    usage
fi

if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: Input BAM not found: $INPUT_BAM"
    exit 1
fi

echo "=========================================="
echo "Downsampling BAM"
echo "=========================================="
echo "Input:  $INPUT_BAM"
echo "Output: $OUTPUT_BAM"
echo "Target: $TARGET_READS reads"
echo "Seed:   $SEED"
echo ""

# Count total reads in input BAM
echo "Counting reads in input BAM..."
TOTAL_READS=$(samtools view -c -@ "$THREADS" "$INPUT_BAM")
echo "Total reads: $TOTAL_READS"

if [ "$TOTAL_READS" -le "$TARGET_READS" ]; then
    echo "Input has fewer reads than target. Copying without downsampling."
    cp "$INPUT_BAM" "$OUTPUT_BAM"
    if [ -f "${INPUT_BAM}.bai" ]; then
        cp "${INPUT_BAM}.bai" "${OUTPUT_BAM}.bai"
    else
        samtools index -@ "$THREADS" "$OUTPUT_BAM"
    fi
    echo "Done."
    exit 0
fi

# Calculate fraction to keep
# samtools -s expects SEED.FRACTION format
# We use bc for floating point division
FRACTION=$(echo "scale=6; $TARGET_READS / $TOTAL_READS" | bc)
echo "Downsampling fraction: $FRACTION"

# Downsample using samtools view -s
echo "Downsampling..."
samtools view -@ "$THREADS" -s "${SEED}.${FRACTION#*.}" -b "$INPUT_BAM" > "$OUTPUT_BAM"

# Index the output
echo "Indexing output BAM..."
samtools index -@ "$THREADS" "$OUTPUT_BAM"

# Verify read count
FINAL_READS=$(samtools view -c -@ "$THREADS" "$OUTPUT_BAM")
echo ""
echo "=========================================="
echo "Downsampling complete"
echo "=========================================="
echo "Original reads: $TOTAL_READS"
echo "Target reads:   $TARGET_READS"
echo "Final reads:    $FINAL_READS"
echo "Output: $OUTPUT_BAM"
