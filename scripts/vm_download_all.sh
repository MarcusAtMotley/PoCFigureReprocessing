#!/bin/bash
# Download ALL files needed for P1, P2, P3, P5 pipelines
# (P4 is already complete)

set -e

echo "=========================================="
echo "Downloading all pipeline inputs"
echo "=========================================="

# ==========================================
# REFERENCES
# ==========================================
echo ""
echo "=== Downloading Reference Files ==="
mkdir -p /data/references
cd /data/references

if [ ! -f "GRCh38_full_analysis_set_plus_decoy_hla.fa" ]; then
    echo "Downloading genome FASTA..."
    aws s3 cp s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa .
    aws s3 cp s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai .
else
    echo "Genome FASTA already exists, skipping..."
fi

if [ ! -d "biscuit_index" ] || [ -z "$(ls -A biscuit_index 2>/dev/null)" ]; then
    echo "Downloading Biscuit index..."
    mkdir -p biscuit_index
    aws s3 cp s3://motleybio/Resources/biscuit_reference_genome/ ./biscuit_index/ --recursive
else
    echo "Biscuit index already exists, skipping..."
fi

# ==========================================
# DNA FASTQs (for P1, P2, P3)
# ==========================================
echo ""
echo "=== Downloading DNA FASTQs ==="
mkdir -p /data/fastq
cd /data/fastq

echo "--- SingleAnalyte WGS (P1) ---"
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R1.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R2.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R1.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R2.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R1.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R2.fastq.gz . --no-progress

echo "--- SingleAnalyte WGEM (P2) ---"
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R1.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R2.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R1.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R2.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R1.fastq.gz . --no-progress
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R2.fastq.gz . --no-progress

echo "--- TrinitySeq (premerged, P1|P2) ---"
aws s3 cp s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/ . --recursive --no-progress

echo "--- HairyTNA (P1|P2) ---"
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21S_3A2_bsTNA-HP_S7_unbarcoded.cutadapt.fastq . --no-progress
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21T_3A2_bsTNA-HP_S8_unbarcoded.cutadapt.fastq . --no-progress
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21W_3A2_bsDNA-HP_S11_unbarcoded.cutadapt.fastq . --no-progress
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21X_3A2_bsDNA-HP_S12_unbarcoded.cutadapt.fastq . --no-progress

# ==========================================
# RNA BAMs (for P5) - from completed P4
# ==========================================
echo ""
echo "=== Downloading RNA BAMs from P4 ==="
mkdir -p /data/rna_bams
cd /data/rna_bams

# Download all P4 BAM outputs
for sample in CoB_02W_1A3_1RNA CoM_02V_1A3_1RNA HT29_02T_1A3_1RNA \
              CoB_08R_3A2_TNA-mRT-EM CoM_08S_3A2_TNA-mRT-EM \
              HT29_21S_3A2_bsTNA-HP HT29_21T_3A2_bsTNA-HP \
              HT29_21U_3A2_bsRNA-HP HT29_21V_3A2_bsRNA-HP; do
    echo "Downloading ${sample}..."
    aws s3 cp "s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p4_rna_counts/${sample}/${sample}.Aligned.sortedByCoord.out.bam" . --no-progress
done

echo ""
echo "=========================================="
echo "Download Complete!"
echo "=========================================="
echo ""
echo "Disk usage:"
du -sh /data/references /data/fastq /data/rna_bams
df -h /data
