#!/bin/bash
# Monitor P3 exon-filtered Batch jobs and post Slack updates
# Posts when: jobs start running, jobs complete, jobs fail, or all done
# Checks every 5 minutes
set -euo pipefail

SLACK_WEBHOOK="${SLACK_WEBHOOK_URL:?Set SLACK_WEBHOOK_URL}"
REGION="us-east-2"
PREFIX="p3ef-"
LOG="/tmp/p3ef_monitor.log"
STATE_FILE="/tmp/p3ef_monitor_state.txt"

TOTAL_JOBS=16

post_slack() {
    local msg="$1"
    echo "$msg"
    echo "[$(date -u '+%H:%M UTC')] $msg" >> "$LOG"
    local payload
    payload=$(python3 -c "import json; print(json.dumps({'text': '''$msg'''}))")
    curl -s -X POST -H 'Content-type: application/json' --data "$payload" "$SLACK_WEBHOOK" > /dev/null 2>&1
}

get_jobs_by_status() {
    local status="$1"
    aws batch list-jobs --job-queue C4_QUEUE --job-status "$status" --region "$REGION" \
        --query "jobSummaryList[?starts_with(jobName, \`$PREFIX\`)].jobName" --output text 2>/dev/null || echo ""
}

count_jobs() {
    local jobs="$1"
    if [ -z "$jobs" ] || [ "$jobs" = "None" ]; then
        echo 0
    else
        echo "$jobs" | wc -w
    fi
}

# Initialize state
touch "$STATE_FILE"
PREV_SUCCEEDED=""
PREV_FAILED=""
PREV_RUNNING=""

post_slack "🔬 *P3 CNV Exon-Filtered Monitor Started*
Tracking $TOTAL_JOBS jobs (prefix: $PREFIX)
Checking every 5 minutes — will notify on completions and failures."

while true; do
    RUNNABLE=$(get_jobs_by_status RUNNABLE)
    STARTING=$(get_jobs_by_status STARTING)
    RUNNING=$(get_jobs_by_status RUNNING)
    SUCCEEDED=$(get_jobs_by_status SUCCEEDED)
    FAILED=$(get_jobs_by_status FAILED)

    N_RUNNABLE=$(count_jobs "$RUNNABLE")
    N_STARTING=$(count_jobs "$STARTING")
    N_RUNNING=$(count_jobs "$RUNNING")
    N_SUCCEEDED=$(count_jobs "$SUCCEEDED")
    N_FAILED=$(count_jobs "$FAILED")
    N_DONE=$((N_SUCCEEDED + N_FAILED))

    # Detect newly completed jobs
    for job in $SUCCEEDED; do
        if ! echo "$PREV_SUCCEEDED" | grep -qw "$job"; then
            # Get filter stats from S3
            SAMPLE="${job#$PREFIX}"
            STATS=$(aws s3 cp "s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p3_cnv_exon_filtered/${SAMPLE}/${SAMPLE}.filter_stats.txt" - 2>/dev/null || echo "")
            PCT=$(echo "$STATS" | grep pct_exonic_removed | cut -d= -f2 || echo "?")
            SEGS=$(echo "$STATS" | grep total_segments | cut -d= -f2 || echo "?")
            post_slack "✅ *$job* completed — ${PCT}% exonic reads removed, ${SEGS} segments
($N_SUCCEEDED/$TOTAL_JOBS succeeded)"
        fi
    done

    # Detect newly failed jobs
    for job in $FAILED; do
        if ! echo "$PREV_FAILED" | grep -qw "$job"; then
            # Try to get failure reason
            JOB_ID=$(aws batch list-jobs --job-queue C4_QUEUE --job-status FAILED --region "$REGION" \
                --query "jobSummaryList[?jobName=='$job'].jobId" --output text 2>/dev/null)
            REASON=$(aws batch describe-jobs --jobs "$JOB_ID" --region "$REGION" \
                --query 'jobs[0].container.reason' --output text 2>/dev/null || echo "unknown")
            post_slack "❌ *$job* FAILED — ${REASON:0:100}
($N_FAILED failed, $N_SUCCEEDED/$TOTAL_JOBS succeeded)"
        fi
    done

    # Detect jobs starting to run (first time we see RUNNING jobs)
    if [ -n "$RUNNING" ] && [ -z "$PREV_RUNNING" ]; then
        post_slack "🚀 Jobs now running: $N_RUNNING active, $N_RUNNABLE queued"
    fi

    # Check if all done
    if [ "$N_DONE" -ge "$TOTAL_JOBS" ]; then
        post_slack "🏁 *All $TOTAL_JOBS P3 exon-filtered jobs complete!*
✅ Succeeded: $N_SUCCEEDED
❌ Failed: $N_FAILED
Results: \`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p3_cnv_exon_filtered/\`"
        break
    fi

    # Save state
    PREV_SUCCEEDED="$SUCCEEDED"
    PREV_FAILED="$FAILED"
    PREV_RUNNING="$RUNNING"

    sleep 300  # 5 minutes
done

echo "Monitor finished at $(date -u)"
