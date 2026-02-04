#!/bin/bash
# Upload all results to S3

set -e

S3_BASE="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"

echo "=========================================="
echo "Uploading Results to S3"
echo "=========================================="

if [ -d "/data/results/p1_p2" ]; then
    echo "Uploading P1+P2 results..."
    aws s3 sync /data/results/p1_p2/ ${S3_BASE}/p1_p2_results/ --no-progress
fi

if [ -d "/data/results/p3" ]; then
    echo "Uploading P3 results..."
    aws s3 sync /data/results/p3/ ${S3_BASE}/p3_cnv_results/ --no-progress
fi

if [ -d "/data/results/p5" ]; then
    echo "Uploading P5 results..."
    aws s3 sync /data/results/p5/ ${S3_BASE}/p5_rna_snp_results/ --no-progress
fi

echo ""
echo "=========================================="
echo "Upload Complete!"
echo "=========================================="
echo ""
echo "Results available at:"
echo "  P1+P2: ${S3_BASE}/p1_p2_results/"
echo "  P3:    ${S3_BASE}/p3_cnv_results/"
echo "  P5:    ${S3_BASE}/p5_rna_snp_results/"
