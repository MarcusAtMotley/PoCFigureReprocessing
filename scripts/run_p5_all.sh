#!/bin/bash
# Run P5 RNA SNP on ALL samples (15 total)
# - 12 TrinitySeq samples (Revelio applied)
# - 3 SingleAnalyte samples (downsampled, no Revelio)

set -e

SCRIPT_DIR=$(dirname "$0")

echo "=========================================="
echo "P5 RNA SNP - ALL Samples"
echo "=========================================="
echo ""
echo "Processing 15 samples total:"
echo ""
echo "TrinitySeq (12 samples) - with Revelio:"
echo "  CoB: TNA-mRT-EM, TNA-RT-EM, RNA-mRT-EM, RNA-RT-EM"
echo "  CoM: TNA-mRT-EM, TNA-RT-EM, RNA-mRT-EM, RNA-RT-EM"
echo "  HT29: bsTNA-HP x2, bsRNA-HP x2"
echo ""
echo "SingleAnalyte (3 samples) - downsampled to 30M, no Revelio:"
echo "  CoB_02W_1A3_1RNA, CoM_02V_1A3_1RNA, HT29_02T_1A3_1RNA"
echo ""

# Check dependencies
echo "Checking dependencies..."
for cmd in samtools bcftools bgzip tabix python3 aws bc; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd not found"
        exit 1
    fi
done
echo "All dependencies found."
echo ""

# Run the P5 pipeline with all samples
exec "$SCRIPT_DIR/simple_p5_rna_snp.sh"
