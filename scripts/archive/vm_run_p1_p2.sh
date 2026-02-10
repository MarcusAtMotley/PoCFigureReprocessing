#!/bin/bash
# Run P1+P2 Pipeline on VM
# Usage: ./vm_run_p1_p2.sh [all|sample_name]

set -e

SAMPLE=${1:-all}

cd /data/PoCFigureReprocessing

# Common Nextflow options - NO scatter-gather, use local Docker
NF_OPTS="-profile docker \
    --genome_fasta /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --genome_fai /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
    --biscuit_index /data/references/biscuit_index \
    --biscuit_chunk_size 0 \
    --revelio_chunk_size 0 \
    -resume"

if [ "$SAMPLE" = "all" ]; then
    echo "=== Running ALL 16 samples ==="
    nextflow run main.nf \
        $NF_OPTS \
        --input /data/PoCFigureReprocessing/samplesheets/p1_p2_local.csv \
        --outdir /data/results
else
    echo "=== Running single sample: $SAMPLE ==="

    # Create single-sample samplesheet
    head -1 /data/PoCFigureReprocessing/samplesheets/p1_p2_local.csv > /tmp/${SAMPLE}.csv
    grep "^${SAMPLE}," /data/PoCFigureReprocessing/samplesheets/p1_p2_local.csv >> /tmp/${SAMPLE}.csv

    nextflow run main.nf \
        $NF_OPTS \
        --input /tmp/${SAMPLE}.csv \
        --outdir /data/results/${SAMPLE}
fi

echo ""
echo "=== Pipeline Complete ==="
echo "Results in: /data/results/"
echo ""
echo "To upload results to S3:"
echo "  aws s3 sync /data/results/ s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p1_p2_results/"
