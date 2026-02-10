#!/bin/bash
# Download all FASTQ files to VM

set -e

echo "=== Downloading FASTQ Files ==="

mkdir -p /data/fastq
cd /data/fastq

echo "--- SingleAnalyte WGS (P1) ---"
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R2.fastq.gz .

echo "--- SingleAnalyte WGEM (P2) ---"
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R2.fastq.gz .

echo "--- TrinitySeq (premerged, P1|P2) ---"
aws s3 cp s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/ . --recursive

echo "--- HairyTNA (P1|P2) ---"
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21S_3A2_bsTNA-HP_S7_unbarcoded.cutadapt.fastq .
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21T_3A2_bsTNA-HP_S8_unbarcoded.cutadapt.fastq .
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21W_3A2_bsDNA-HP_S11_unbarcoded.cutadapt.fastq .
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21X_3A2_bsDNA-HP_S12_unbarcoded.cutadapt.fastq .

echo ""
echo "=== FASTQs Downloaded ==="
ls -lh /data/fastq/
du -sh /data/fastq/
