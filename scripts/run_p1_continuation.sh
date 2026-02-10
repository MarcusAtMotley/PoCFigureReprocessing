#!/bin/bash
# P1 DNA SNP Continuation Script
# Monitors batch 1, then automatically starts batch 2
set -e

SCRIPT_DIR=$(dirname "$0")
RESULTS=/data/results/p1_dna_snp

echo "=========================================="
echo "P1 DNA SNP - Continuation Monitor"
echo "=========================================="
echo "Started monitoring: $(date)"
echo ""

# Batch 1 samples to monitor
BATCH1_SAMPLES=(
    "CoB_08L_3A2_DNA-EM"
    "CoM_08M_3A2_DNA-EM"
    "CoB_08R_3A2_TNA-mRT-EM"
)

# Function to check if a sample is complete (has VCF output)
check_sample_complete() {
    local SAMPLE=$1
    if [ -f "$RESULTS/${SAMPLE}/${SAMPLE}.bcftools.vcf" ]; then
        return 0
    fi
    return 1
}

# Wait for all batch 1 samples to complete
echo "Waiting for Batch 1 to complete..."
while true; do
    COMPLETE=0
    for SAMPLE in "${BATCH1_SAMPLES[@]}"; do
        if check_sample_complete "$SAMPLE"; then
            ((COMPLETE++))
        fi
    done

    echo "[$(date '+%H:%M')] Batch 1 progress: $COMPLETE/${#BATCH1_SAMPLES[@]} samples complete"

    if [ $COMPLETE -eq ${#BATCH1_SAMPLES[@]} ]; then
        echo ""
        echo "=========================================="
        echo "Batch 1 Complete! Starting Batch 2..."
        echo "=========================================="
        break
    fi

    # Check every 15 minutes
    sleep 900
done

# Start batch 2
echo ""
echo "Launching Batch 2..."
"$SCRIPT_DIR/run_p1_batch2.sh"

echo ""
echo "=========================================="
echo "All P1 DNA SNP batches complete!"
echo "=========================================="
echo "Finished: $(date)"
