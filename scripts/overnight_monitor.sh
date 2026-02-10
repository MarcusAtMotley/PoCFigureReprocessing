#!/bin/bash
# Overnight monitoring script for P1 DNA SNP jobs
# Checks Batch job status + local CoM_08M progress every 30 minutes
# Writes to /data/overnight_monitor.log

LOG="/data/overnight_monitor.log"
S3_OUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p1_dna_snp"

echo "========================================" >> "$LOG"
echo "Overnight Monitor Started: $(date)" >> "$LOG"
echo "========================================" >> "$LOG"

while true; do
    echo "" >> "$LOG"
    echo "--- Check at $(date) ---" >> "$LOG"

    # Check Batch jobs
    echo "" >> "$LOG"
    echo "[Batch Jobs]" >> "$LOG"
    for STATUS in SUCCEEDED RUNNING RUNNABLE FAILED; do
        JOBS=$(aws batch list-jobs --job-queue C4_QUEUE --job-status $STATUS \
            --query 'jobSummaryList[*].jobName' --output text 2>/dev/null)
        if [ -n "$JOBS" ]; then
            echo "  $STATUS: $JOBS" >> "$LOG"
        fi
    done

    # Check local CoM_08M progress
    echo "" >> "$LOG"
    echo "[Local: CoM_08M_3A2_DNA-EM]" >> "$LOG"
    if ls /data/dna_work/CoM_08M_3A2_DNA-EM/*.bcftools.vcf 2>/dev/null; then
        VARIANTS=$(grep -vc '^#' /data/results/p1_dna_snp/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf 2>/dev/null || echo "?")
        echo "  STATUS: COMPLETE ($VARIANTS variants)" >> "$LOG"
    elif ls /data/dna_work/CoM_08M_3A2_DNA-EM/*.revelio.sorted.bam 2>/dev/null | grep -qv tmp; then
        echo "  STATUS: Revelio done, BCFtools running" >> "$LOG"
    elif ls /data/dna_work/CoM_08M_3A2_DNA-EM/*.calmd.bam 2>/dev/null; then
        echo "  STATUS: Calmd done, Revelio running" >> "$LOG"
    elif ls /data/dna_work/CoM_08M_3A2_DNA-EM/*.markdup.bam 2>/dev/null; then
        echo "  STATUS: Markdup done, Calmd running" >> "$LOG"
    elif ls /data/dna_work/CoM_08M_3A2_DNA-EM/*.sorted.bam 2>/dev/null | grep -qv tmp; then
        echo "  STATUS: Alignment done, Markdup running" >> "$LOG"
    else
        TMPS=$(ls /data/dna_work/CoM_08M_3A2_DNA-EM/*.sorted.bam.tmp.*.bam 2>/dev/null | wc -l)
        echo "  STATUS: Still aligning ($TMPS sort temp files)" >> "$LOG"
    fi

    # If CoM_08M finishes, auto-upload to S3
    if [ -f "/data/results/p1_dna_snp/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf" ] && \
       [ ! -f "/data/dna_work/CoM_08M_uploaded" ]; then
        echo "  Uploading CoM_08M results to S3..." >> "$LOG"
        aws s3 cp /data/results/p1_dna_snp/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf "$S3_OUT/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf" --quiet
        aws s3 cp /data/results/p1_dna_snp/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf.gz "$S3_OUT/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf.gz" --quiet
        aws s3 cp /data/results/p1_dna_snp/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf.gz.tbi "$S3_OUT/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf.gz.tbi" --quiet
        aws s3 cp /data/dna_work/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.revelio.sorted.bam "$S3_OUT/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.revelio.sorted.bam" --quiet
        aws s3 cp /data/dna_work/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.revelio.sorted.bam.bai "$S3_OUT/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.revelio.sorted.bam.bai" --quiet
        aws s3 cp /data/dna_work/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.markdup_metrics.txt "$S3_OUT/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.markdup_metrics.txt" --quiet
        touch /data/dna_work/CoM_08M_uploaded
        echo "  Upload complete!" >> "$LOG"
    fi

    # Check if everything is done
    BATCH_RUNNING=$(aws batch list-jobs --job-queue C4_QUEUE --job-status RUNNING --query 'length(jobSummaryList)' --output text 2>/dev/null)
    BATCH_RUNNABLE=$(aws batch list-jobs --job-queue C4_QUEUE --job-status RUNNABLE --query 'length(jobSummaryList)' --output text 2>/dev/null)
    LOCAL_DONE=false
    [ -f "/data/results/p1_dna_snp/CoM_08M_3A2_DNA-EM/CoM_08M_3A2_DNA-EM.bcftools.vcf" ] && LOCAL_DONE=true

    if [ "$BATCH_RUNNING" = "0" ] && [ "$BATCH_RUNNABLE" = "0" ] && [ "$LOCAL_DONE" = "true" ]; then
        echo "" >> "$LOG"
        echo "========================================" >> "$LOG"
        echo "ALL P1 DNA SNP PROCESSING COMPLETE!" >> "$LOG"
        echo "$(date)" >> "$LOG"
        echo "========================================" >> "$LOG"

        # Final summary
        echo "" >> "$LOG"
        echo "[Final Summary]" >> "$LOG"
        SUCCEEDED=$(aws batch list-jobs --job-queue C4_QUEUE --job-status SUCCEEDED --query 'jobSummaryList[*].jobName' --output text 2>/dev/null)
        FAILED=$(aws batch list-jobs --job-queue C4_QUEUE --job-status FAILED --query 'jobSummaryList[*].jobName' --output text 2>/dev/null)
        echo "  Batch SUCCEEDED: $SUCCEEDED" >> "$LOG"
        [ -n "$FAILED" ] && echo "  Batch FAILED: $FAILED" >> "$LOG"
        echo "  Local: CoB_08L DONE, CoB_08R DONE, CoM_08M DONE" >> "$LOG"

        # List all S3 results
        echo "" >> "$LOG"
        echo "[S3 Results]" >> "$LOG"
        aws s3 ls "$S3_OUT/" --recursive 2>/dev/null | grep -E "\.vcf$|\.vcf\.gz$" >> "$LOG"
        break
    fi

    sleep 1800  # 30 minutes
done
