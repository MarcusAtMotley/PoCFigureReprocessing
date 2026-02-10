#!/bin/bash
# Run P5 RNA SNP on all TrinitySeq samples (12 total)
# This uses existing P4 BAMs from S3
#
# Execution time estimate: ~2-4 hours for all 12 samples with BCFtools
# (compared to 16+ hours per sample with LoFreq)

set -e

SCRIPT_DIR=$(dirname "$0")

echo "=========================================="
echo "P5 RNA SNP - TrinitySeq Samples"
echo "=========================================="
echo ""
echo "This script processes 12 TrinitySeq RNA samples through:"
echo "  1. Download P4 BAM from S3"
echo "  2. samtools calmd (add MD tags)"
echo "  3. Revelio (mask bisulfite/EM artifacts)"
echo "  4. BCFtools (variant calling)"
echo ""
echo "Samples:"
echo "  CoB: 4 mTNA samples (mRT, RT, mRNA-mRT, mRNA-RT)"
echo "  CoM: 4 mTNA samples (mRT, RT, mRNA-mRT, mRNA-RT)"
echo "  HT29: 4 HairyTNA samples (bsTNA x2, bsRNA x2)"
echo ""

# Check dependencies
echo "Checking dependencies..."
for cmd in samtools bcftools bgzip tabix python3 aws; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd not found"
        exit 1
    fi
done
echo "All dependencies found."
echo ""

# Run the P5 pipeline with TrinitySeq samples only
exec "$SCRIPT_DIR/simple_p5_rna_snp.sh" --trinityseq-only
