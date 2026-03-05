#!/bin/bash
# Submit P3 CNV exon-filtered re-run to AWS Batch (all 16 samples in parallel)
#
# Usage: submit_p3_exon_filtered.sh [--dry-run]
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
P3_SAMPLESHEET="${SCRIPT_DIR}/../samplesheets/p3_batch.csv"
JOB_QUEUE="C4_QUEUE"
JOB_DEF="poc-dna-pipeline"

DRY_RUN=false
[ "${1:-}" = "--dry-run" ] && DRY_RUN=true

echo "=========================================="
echo "Submitting P3 CNV Exon-Filtered Batch Jobs"
echo "=========================================="
echo "Queue: $JOB_QUEUE"
echo "Samplesheet: $P3_SAMPLESHEET"
echo "Output: s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p3_cnv_exon_filtered/"
[ "$DRY_RUN" = true ] && echo "[DRY RUN MODE]"
echo ""

SUBMITTED=0
JOB_IDS=""

tail -n +2 "$P3_SAMPLESHEET" | while IFS=',' read -r SAMPLE BAM_PATH BAM_TYPE; do
    [ -z "$SAMPLE" ] && continue

    # WGS sorted BAMs need markdup (more resources + time)
    if [ "$BAM_TYPE" = "sorted" ]; then
        VCPUS="16"
        MEMORY="60000"
    else
        VCPUS="16"
        MEMORY="30000"
    fi

    echo "Submitting: $SAMPLE ($BAM_TYPE BAM)"
    echo "  BAM: ${BAM_PATH:0:80}..."
    echo "  Resources: ${VCPUS} vCPU, ${MEMORY}MB RAM"

    CMD_ARGS="[\"/opt/batch_p3_cnv_exon_filtered.sh\", \"$SAMPLE\", \"$BAM_PATH\", \"$BAM_TYPE\"]"

    if [ "$DRY_RUN" = true ]; then
        echo "  [DRY RUN] Would submit: $CMD_ARGS"
    else
        JOB_ID=$(aws batch submit-job \
            --job-name "p3ef-${SAMPLE}" \
            --job-queue "$JOB_QUEUE" \
            --job-definition "$JOB_DEF" \
            --container-overrides "{
                \"command\": $CMD_ARGS,
                \"resourceRequirements\": [
                    {\"type\": \"VCPU\", \"value\": \"$VCPUS\"},
                    {\"type\": \"MEMORY\", \"value\": \"$MEMORY\"}
                ]
            }" \
            --query 'jobId' --output text 2>&1)
        echo "  Job ID: $JOB_ID"
        JOB_IDS="$JOB_IDS $JOB_ID"
    fi

    echo ""
    SUBMITTED=$((SUBMITTED+1))
done

echo "=========================================="
echo "Submitted $SUBMITTED jobs"
echo "=========================================="
echo ""
echo "Monitor with:"
echo "  aws batch list-jobs --job-queue $JOB_QUEUE --job-status RUNNING --query 'jobSummaryList[?starts_with(jobName,\`p3ef-\`)].{name:jobName,status:status,id:jobId}' --output table"
echo "  aws batch list-jobs --job-queue $JOB_QUEUE --job-status SUCCEEDED --query 'jobSummaryList[?starts_with(jobName,\`p3ef-\`)].{name:jobName,status:status}' --output table"
