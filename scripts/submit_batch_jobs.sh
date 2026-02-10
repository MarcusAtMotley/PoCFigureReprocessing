#!/bin/bash
# Submit P1 DNA SNP jobs to AWS Batch using the samplesheet
# Each sample runs on its own on-demand instance
set -e

SAMPLESHEET="${1:-$(dirname $0)/../samplesheets/p1_dna_snp.csv}"
JOB_QUEUE="C4_QUEUE"
JOB_DEF="poc-dna-pipeline"

# Samples to skip (already processing locally)
SKIP_SAMPLES="CoB_08L_3A2_DNA-EM CoM_08M_3A2_DNA-EM CoB_08R_3A2_TNA-mRT-EM"

echo "=========================================="
echo "Submitting P1 DNA SNP Batch Jobs"
echo "=========================================="
echo "Queue: $JOB_QUEUE"
echo "Samplesheet: $SAMPLESHEET"
echo "Skipping: $SKIP_SAMPLES"
echo ""

SUBMITTED=0
SKIPPED=0

# Read samplesheet (skip header)
tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r SAMPLE FQ1 FQ2 SINGLE_END PIPELINE CELL_LINE ASSAY_CATEGORY; do
    # Skip empty lines
    [ -z "$SAMPLE" ] && continue

    # Skip samples already running locally
    if echo "$SKIP_SAMPLES" | grep -qw "$SAMPLE"; then
        echo "SKIP: $SAMPLE (running locally)"
        SKIPPED=$((SKIPPED+1))
        continue
    fi

    # Determine resource requirements based on assay type
    # SingleAnalyte WGS has huge files, HairyTNA is small single-end
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

    echo "Submitting: $SAMPLE ($ASSAY_CATEGORY, ${SINGLE_END})"
    echo "  R1: ${FQ1:0:80}..."
    echo "  R2: ${FQ2:0:80}..."
    echo "  Resources: ${VCPUS} vCPU, ${MEMORY}MB RAM"

    # Build command args
    CMD_ARGS="[\"/opt/batch_dna_pipeline.sh\", \"$SAMPLE\", \"$SINGLE_END\", \"$FQ1\""
    if [ -n "$FQ2" ]; then
        CMD_ARGS="$CMD_ARGS, \"$FQ2\""
    fi
    CMD_ARGS="$CMD_ARGS]"

    JOB_ID=$(aws batch submit-job \
        --job-name "p1-${SAMPLE}" \
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
    echo ""
    SUBMITTED=$((SUBMITTED+1))
done

echo "=========================================="
echo "Done! Submitted: $SUBMITTED, Skipped: $SKIPPED"
echo "=========================================="
echo ""
echo "Monitor with:"
echo "  aws batch list-jobs --job-queue $JOB_QUEUE --job-status RUNNING"
echo "  aws batch describe-jobs --jobs <JOB_ID>"
