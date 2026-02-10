#!/bin/bash
# Run P1 DNA SNP Pipeline
# Uses simple_dna_pipeline.sh with parallel Revelio + BCFtools

set -e

SCRIPT_DIR=$(dirname "$0")
source "$SCRIPT_DIR/simple_dna_pipeline.sh"

echo "=========================================="
echo "P1 DNA SNP Pipeline"
echo "=========================================="
echo ""

# Check dependencies
echo "Checking dependencies..."
for cmd in samtools bcftools bgzip tabix python3 bwa-mem2; do
    if ! command -v $cmd &> /dev/null; then
        echo "WARNING: $cmd not found (may use docker instead)"
    fi
done
echo ""

# Process available samples
echo "Processing DNA samples..."
echo ""

# CoB DNA-EM (local FASTQ available)
if [ -f "/data/fastq/CoB_08L_3A2_DNA-EM_S1_merged_R1.fastq" ]; then
    echo "=== CoB_08L_3A2_DNA-EM ==="
    process_dna_sample "CoB_08L_3A2_DNA-EM" \
        "/data/fastq/CoB_08L_3A2_DNA-EM_S1_merged_R1.fastq" \
        "/data/fastq/CoB_08L_3A2_DNA-EM_S1_merged_R2.fastq" \
        "false"
fi

echo ""
echo "=========================================="
echo "P1 DNA SNP Complete!"
echo "=========================================="
echo "Results in: $RESULTS/p1_dna_snp/"
ls -la "$RESULTS/p1_dna_snp/" 2>/dev/null || echo "No results yet"
