#!/bin/bash
# 3-hour status check — queries all running jobs and posts summary to Slack
set -euo pipefail

SLACK_WEBHOOK="${SLACK_WEBHOOK_URL:?Set SLACK_WEBHOOK_URL environment variable}"
REGION="us-east-2"

echo "=== Status check starting at $(date -u) ==="

# Get running jobs
RUNNING=$(aws batch list-jobs --job-queue C4_QUEUE --job-status RUNNING --region $REGION 2>/dev/null)
RUNNING_COUNT=$(echo "$RUNNING" | python3 -c "import json,sys; data=json.load(sys.stdin); print(len(data['jobSummaryList']))")

# Get recently succeeded jobs (last 24h)
SUCCEEDED=$(aws batch list-jobs --job-queue C4_QUEUE --job-status SUCCEEDED --region $REGION 2>/dev/null)
FAILED=$(aws batch list-jobs --job-queue C4_QUEUE --job-status FAILED --region $REGION 2>/dev/null)

# Build running job details
RUNNING_DETAILS=$(echo "$RUNNING" | python3 -c "
import json, sys
from datetime import datetime, timezone
data = json.load(sys.stdin)
now = datetime.now(timezone.utc).timestamp()
lines = []
for j in data['jobSummaryList']:
    hours = (now - j['startedAt']/1000) / 3600
    lines.append(f\"  • {j['jobName']} — {hours:.1f}h running\")
print('\n'.join(lines) if lines else '  None')
")

# Check for newly completed jobs (stopped in last 3h)
NEW_COMPLETED=$(echo "$SUCCEEDED" | python3 -c "
import json, sys
from datetime import datetime, timezone
data = json.load(sys.stdin)
now = datetime.now(timezone.utc).timestamp()
lines = []
for j in data['jobSummaryList']:
    stopped = j.get('stoppedAt', 0) / 1000
    if now - stopped < 10800:  # 3 hours
        hours = (stopped - j['createdAt']/1000) / 3600
        lines.append(f\"  ✅ {j['jobName']} — {hours:.1f}h total\")
print('\n'.join(lines) if lines else '  None')
")

# Check for newly failed jobs (stopped in last 3h)
NEW_FAILED=$(echo "$FAILED" | python3 -c "
import json, sys
from datetime import datetime, timezone
data = json.load(sys.stdin)
now = datetime.now(timezone.utc).timestamp()
lines = []
for j in data['jobSummaryList']:
    stopped = j.get('stoppedAt', 0) / 1000
    if now - stopped < 10800:  # 3 hours
        reason = j.get('container',{}).get('reason','unknown')
        lines.append(f\"  ❌ {j['jobName']} — {reason[:60]}\")
print('\n'.join(lines) if lines else '  None')
")

# Get latest log line for each running job
LOG_DETAILS=$(echo "$RUNNING" | python3 -c "
import json, sys, subprocess
data = json.load(sys.stdin)
for j in data['jobSummaryList']:
    jid = j['jobId']
    try:
        desc = json.loads(subprocess.check_output(
            ['aws', 'batch', 'describe-jobs', '--jobs', jid, '--region', 'us-east-2'],
            stderr=subprocess.DEVNULL).decode())
        ls = desc['jobs'][0]['container']['logStreamName']
        logs = json.loads(subprocess.check_output(
            ['aws', 'logs', 'get-log-events', '--log-group-name', '/aws/batch/job',
             '--log-stream-name', ls, '--limit', '1', '--no-start-from-head',
             '--region', 'us-east-2'], stderr=subprocess.DEVNULL).decode())
        if logs['events']:
            msg = logs['events'][0]['message'][:80]
            print(f\"  📋 {j['jobName']}: {msg}\")
    except:
        print(f\"  📋 {j['jobName']}: (could not fetch logs)\")
")

TIMESTAMP=$(date -u '+%Y-%m-%d %H:%M UTC')

# Compose Slack message
MESSAGE="🧬 *Pipeline Status Check — ${TIMESTAMP}*

*Running (${RUNNING_COUNT} jobs):*
${RUNNING_DETAILS}

*Latest log per job:*
${LOG_DETAILS}

*Completed in last 3h:*
${NEW_COMPLETED}

*Failed in last 3h:*
${NEW_FAILED}"

echo "$MESSAGE"

# Post to Slack
PAYLOAD=$(python3 -c "
import json
msg = '''${MESSAGE}'''
print(json.dumps({'text': msg}))
")

curl -s -X POST -H 'Content-type: application/json' --data "$PAYLOAD" "$SLACK_WEBHOOK"
echo ""
echo "=== Slack notification sent ==="
