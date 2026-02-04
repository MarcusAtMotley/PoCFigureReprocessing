#!/bin/bash
# Download reference files to VM

set -e

echo "=== Downloading Reference Files ==="

mkdir -p /data/references
cd /data/references

echo "Downloading genome FASTA..."
aws s3 cp s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa .
aws s3 cp s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai .

echo "Downloading Biscuit index..."
aws s3 cp s3://motleybio/Resources/biscuit_reference_genome/ ./biscuit_index/ --recursive

echo ""
echo "=== References Downloaded ==="
ls -lh /data/references/
