#!/bin/bash
# Submit P5 RNA SNP re-runs for 6 samples missing .bcftools.vcf.gz
# These samples only have old .vcf.gz from the early Seqera/LoFreq run
#
# Usage: submit_p5_rerun.sh [--dry-run]
set -e

JOB_QUEUE="C4_QUEUE"
JOB_DEF="poc-dna-pipeline"
REGION="us-east-2"
S3_BASE="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p4_rna_counts"

DRY_RUN=false
[ "${1:-}" = "--dry-run" ] && DRY_RUN=true

echo "=========================================="
echo "Submitting P5 RNA SNP Re-runs (bcftools)"
echo "=========================================="
[ "$DRY_RUN" = true ] && echo "[DRY RUN MODE]"
echo ""

# 6 samples that only have old .vcf.gz naming
SAMPLES=(
    "CoB_08X_3A2_TNA-RT-EM"
    "CoB_09D_3A2_RNA-mRT-EM"
    "CoB_09J_3A2_RNA-RT-EM"
    "CoM_08Y_3A2_TNA-RT-EM"
    "CoM_09E_3A2_RNA-mRT-EM"
    "CoM_09K_3A2_RNA-RT-EM"
)

SUBMITTED=0
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_PATH="${S3_BASE}/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    JOB_NAME="p5-rerun-${SAMPLE}"

    echo "Sample: $SAMPLE"
    echo "  BAM: $BAM_PATH"

    if [ "$DRY_RUN" = true ]; then
        echo "  [DRY RUN] Would submit: $JOB_NAME"
    else
        JOB_ID=$(aws batch submit-job \
            --job-name "$JOB_NAME" \
            --job-queue "$JOB_QUEUE" \
            --job-definition "$JOB_DEF" \
            --container-overrides "{
                \"command\": [\"/opt/batch_p5_rna_snp.sh\", \"$SAMPLE\", \"$BAM_PATH\"],
                \"resourceRequirements\": [
                    {\"type\": \"VCPU\", \"value\": \"16\"},
                    {\"type\": \"MEMORY\", \"value\": \"30000\"}
                ]
            }" \
            --region "$REGION" \
            --query 'jobId' --output text)
        echo "  Submitted: $JOB_ID"
        SUBMITTED=$((SUBMITTED + 1))
    fi
    echo ""
done

echo "=========================================="
echo "Submitted $SUBMITTED / ${#SAMPLES[@]} jobs"
echo "=========================================="
