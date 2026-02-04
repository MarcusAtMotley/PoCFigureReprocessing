#!/bin/bash
# Run all pipelines (P1, P2, P3, P5) on VM
# P4 is already complete
#
# Usage:
#   ./vm_run_all_pipelines.sh           # Run all
#   ./vm_run_all_pipelines.sh p1        # Run only P1
#   ./vm_run_all_pipelines.sh p1_p2     # Run P1 and P2
#   ./vm_run_all_pipelines.sh p5        # Run only P5

set -e

PIPELINE=${1:-all}

cd /data/PoCFigureReprocessing

# Common options - NO scatter-gather on VM (local I/O is fast)
BASE_OPTS="-profile docker \
    --genome_fasta /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --genome_fai /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
    --biscuit_chunk_size 0 \
    --revelio_chunk_size 0 \
    -resume"

DNA_OPTS="$BASE_OPTS --biscuit_index /data/references/biscuit_index"

run_p1_p2() {
    echo ""
    echo "=========================================="
    echo "Running P1 (DNA SNP) + P2 (DNA Methylation)"
    echo "=========================================="
    nextflow run main.nf \
        $DNA_OPTS \
        --input /data/PoCFigureReprocessing/samplesheets/p1_p2_local.csv \
        --outdir /data/results/p1_p2
}

run_p3() {
    echo ""
    echo "=========================================="
    echo "Running P3 (CNV)"
    echo "=========================================="
    nextflow run main.nf \
        $DNA_OPTS \
        --input /data/PoCFigureReprocessing/samplesheets/p3_local.csv \
        --outdir /data/results/p3
}

run_p5() {
    echo ""
    echo "=========================================="
    echo "Running P5 (RNA SNP)"
    echo "=========================================="
    nextflow run main.nf \
        $BASE_OPTS \
        --input /data/PoCFigureReprocessing/samplesheets/p5_local.csv \
        --outdir /data/results/p5
}

case $PIPELINE in
    p1_p2)
        run_p1_p2
        ;;
    p3)
        run_p3
        ;;
    p5)
        run_p5
        ;;
    all)
        run_p1_p2
        run_p3
        run_p5
        ;;
    *)
        echo "Usage: $0 [all|p1_p2|p3|p5]"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Pipeline(s) Complete!"
echo "=========================================="
echo ""
echo "Results in: /data/results/"
ls -la /data/results/
echo ""
echo "To upload: ./scripts/vm_upload_results.sh"
