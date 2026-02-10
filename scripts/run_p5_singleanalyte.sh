#!/bin/bash
# Run P5 RNA SNP on SingleAnalyte samples (3 total, downsampled)
# These have much higher depth than TrinitySeq, so we downsample to 30M reads
#
# Note: Revelio is applied to ALL samples for pipeline consistency (comparability)

set -e

SCRIPT_DIR=$(dirname "$0")

echo "=========================================="
echo "P5 RNA SNP - SingleAnalyte Samples"
echo "=========================================="
echo ""
echo "This script processes 3 SingleAnalyte RNA samples through:"
echo "  1. Download P4 BAM from S3"
echo "  2. Downsample to 30M reads (match TrinitySeq depth)"
echo "  3. samtools calmd (add MD tags)"
echo "  4. Revelio (applied to ALL samples for comparability)"
echo "  5. BCFtools (variant calling)"
echo ""
echo "Samples:"
echo "  CoB_02W_1A3_1RNA"
echo "  CoM_02V_1A3_1RNA"
echo "  HT29_02T_1A3_1RNA"
echo ""

# Check dependencies
echo "Checking dependencies..."
for cmd in samtools bcftools bgzip tabix aws bc; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd not found"
        exit 1
    fi
done
echo "All dependencies found."
echo ""

# Run the P5 pipeline with SingleAnalyte samples only
exec "$SCRIPT_DIR/simple_p5_rna_snp.sh" --singleanalyte-only
