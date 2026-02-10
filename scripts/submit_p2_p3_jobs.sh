#!/bin/bash
# Submit P2 DNA Methylation and P3 CNV jobs to AWS Batch
# P2 runs first; P3 uses markdup BAMs from P2 output (or sorted BAMs from P1)
#
# Usage: submit_p2_p3_jobs.sh [p2|p3|both] [--dry-run]
#   p2:       Submit only P2 methylation jobs
#   p3:       Submit only P3 CNV jobs
#   both:     Submit P2 then P3 (default)
#   --dry-run: Print commands without submitting
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
P2_SAMPLESHEET="${SCRIPT_DIR}/../samplesheets/p2_batch.csv"
P3_SAMPLESHEET="${SCRIPT_DIR}/../samplesheets/p3_batch.csv"
JOB_QUEUE="C4_QUEUE"
JOB_DEF="poc-dna-pipeline"

MODE="${1:-both}"
DRY_RUN=false
for arg in "$@"; do
    [ "$arg" = "--dry-run" ] && DRY_RUN=true
done

submit_p2_jobs() {
    echo "=========================================="
    echo "Submitting P2 DNA Methylation Batch Jobs"
    echo "=========================================="
    echo "Queue: $JOB_QUEUE"
    echo "Samplesheet: $P2_SAMPLESHEET"
    echo ""

    SUBMITTED=0

    tail -n +2 "$P2_SAMPLESHEET" | while IFS=',' read -r SAMPLE INPUT_MODE INPUT_1 INPUT_2 SINGLE_END CELL_LINE ASSAY_CATEGORY; do
        [ -z "$SAMPLE" ] && continue

        # Resource sizing: WGEM FASTQ alignment needs more resources
        case "$ASSAY_CATEGORY" in
            single_analyte)
                VCPUS="32"
                MEMORY="60000"
                ;;
            HairyTNA)
                VCPUS="16"
                MEMORY="30000"
                ;;
            *)
                VCPUS="32"
                MEMORY="56000"
                ;;
        esac

        echo "Submitting P2: $SAMPLE ($INPUT_MODE, $ASSAY_CATEGORY)"
        echo "  Input: ${INPUT_1:0:80}..."
        echo "  Resources: ${VCPUS} vCPU, ${MEMORY}MB RAM"

        # Build command args
        CMD_ARGS="[\"/opt/batch_p2_methylation.sh\", \"$SAMPLE\", \"$SINGLE_END\", \"$INPUT_MODE\", \"$INPUT_1\""
        if [ -n "$INPUT_2" ]; then
            CMD_ARGS="$CMD_ARGS, \"$INPUT_2\""
        fi
        CMD_ARGS="$CMD_ARGS]"

        if [ "$DRY_RUN" = true ]; then
            echo "  [DRY RUN] Would submit: $CMD_ARGS"
        else
            JOB_ID=$(aws batch submit-job \
                --job-name "p2-${SAMPLE}" \
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
        fi

        echo ""
        SUBMITTED=$((SUBMITTED+1))
    done

    echo "P2 submission complete ($SUBMITTED samples)"
    echo ""
}

submit_p3_jobs() {
    echo "=========================================="
    echo "Submitting P3 CNV Batch Jobs"
    echo "=========================================="
    echo "Queue: $JOB_QUEUE"
    echo "Samplesheet: $P3_SAMPLESHEET"
    echo ""

    SUBMITTED=0

    tail -n +2 "$P3_SAMPLESHEET" | while IFS=',' read -r SAMPLE BAM_PATH BAM_TYPE; do
        [ -z "$SAMPLE" ] && continue

        # CNVpytor is memory-heavy; WGS samples are larger
        VCPUS="16"
        MEMORY="30000"

        echo "Submitting P3: $SAMPLE ($BAM_TYPE BAM)"
        echo "  BAM: ${BAM_PATH:0:80}..."
        echo "  Resources: ${VCPUS} vCPU, ${MEMORY}MB RAM"

        CMD_ARGS="[\"/opt/batch_p3_cnv.sh\", \"$SAMPLE\", \"$BAM_PATH\", \"$BAM_TYPE\"]"

        if [ "$DRY_RUN" = true ]; then
            echo "  [DRY RUN] Would submit: $CMD_ARGS"
        else
            JOB_ID=$(aws batch submit-job \
                --job-name "p3-${SAMPLE}" \
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
        fi

        echo ""
        SUBMITTED=$((SUBMITTED+1))
    done

    echo "P3 submission complete ($SUBMITTED samples)"
    echo ""
}

case "$MODE" in
    p2)
        submit_p2_jobs
        ;;
    p3)
        submit_p3_jobs
        ;;
    both|*)
        submit_p2_jobs
        submit_p3_jobs
        ;;
esac

echo "=========================================="
echo "Monitor with:"
echo "  aws batch list-jobs --job-queue $JOB_QUEUE --job-status RUNNING"
echo "  aws batch describe-jobs --jobs <JOB_ID>"
echo "=========================================="
