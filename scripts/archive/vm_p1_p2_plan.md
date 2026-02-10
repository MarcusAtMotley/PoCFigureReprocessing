# P1+P2 Pipeline on VM - Execution Plan

Run the full P1 (DNA SNP) and P2 (DNA Methylation) pipeline on a large VM with local NVMe storage to avoid S3 overhead.

## VM Requirements

- **Instance type:** i3.8xlarge or i3.16xlarge (NVMe storage for fast I/O)
- **Storage:** Use instance NVMe (~7-14 TB)
- **Region:** us-east-2 (same as S3 buckets)

## Setup Steps

### 1. Launch VM and Mount NVMe

```bash
# SSH into VM
ssh -i your-key.pem ubuntu@<VM_IP>

# Mount NVMe drives (i3 instances)
sudo mkfs.ext4 /dev/nvme0n1
sudo mkdir /data
sudo mount /dev/nvme0n1 /data
sudo chown ubuntu:ubuntu /data
cd /data
```

### 2. Install Dependencies

```bash
# Install Java (required for Nextflow)
sudo apt update
sudo apt install -y openjdk-17-jdk

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Install Docker
sudo apt install -y docker.io
sudo usermod -aG docker ubuntu
newgrp docker

# Install AWS CLI
sudo apt install -y awscli
aws configure  # Enter credentials
```

### 3. Clone Pipeline

```bash
cd /data
git clone https://github.com/MarcusAtMotley/PoCFigureReprocessing.git
cd PoCFigureReprocessing
```

### 4. Download Reference Files (~30 GB)

```bash
mkdir -p /data/references
cd /data/references

# Genome FASTA + index
aws s3 cp s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa .
aws s3 cp s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai .

# Biscuit index
aws s3 cp s3://motleybio/Resources/biscuit_reference_genome/ ./biscuit_index/ --recursive
```

### 5. Download Input FASTQs

```bash
mkdir -p /data/fastq
cd /data/fastq

# === SingleAnalyte WGS (P1 only) ===
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R2.fastq.gz .

# === SingleAnalyte WGEM (P2 only) ===
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R2.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R1.fastq.gz .
aws s3 cp s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R2.fastq.gz .

# === TrinitySeq (premerged, P1|P2) ===
aws s3 cp s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/ . --recursive

# === HairyTNA (P1|P2) ===
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21S_3A2_bsTNA-HP_S7_unbarcoded.cutadapt.fastq .
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21T_3A2_bsTNA-HP_S8_unbarcoded.cutadapt.fastq .
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21W_3A2_bsDNA-HP_S11_unbarcoded.cutadapt.fastq .
aws s3 cp s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21X_3A2_bsDNA-HP_S12_unbarcoded.cutadapt.fastq .
```

## Running the Pipeline

### Option A: Run All at Once

Create a local samplesheet pointing to `/data/fastq/` files, then run:

```bash
cd /data/PoCFigureReprocessing

nextflow run main.nf \
    -profile docker \
    --input /data/samplesheets/p1_p2_local.csv \
    --outdir /data/results \
    --genome_fasta /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --genome_fai /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
    --biscuit_index /data/references/biscuit_index \
    --biscuit_chunk_size 0 \
    --revelio_chunk_size 0 \
    -resume
```

### Option B: Run Sample by Sample (Safer)

Process one sample at a time to avoid memory issues:

```bash
# Script: run_single_sample.sh
SAMPLE=$1
PIPELINE=$2  # P1, P2, or P1|P2

nextflow run main.nf \
    -profile docker \
    --input /data/samplesheets/${SAMPLE}.csv \
    --outdir /data/results/${SAMPLE} \
    --genome_fasta /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --genome_fai /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
    --biscuit_index /data/references/biscuit_index \
    --biscuit_chunk_size 0 \
    --revelio_chunk_size 0 \
    -resume
```

Run samples sequentially:
```bash
./run_single_sample.sh CoB_02M_1C3_1DNA P1
./run_single_sample.sh CoM_02K_1C3_1DNA P1
# ... etc
```

## Upload Results

```bash
# Upload results to S3
aws s3 sync /data/results/ s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p1_p2_results/
```

## Estimated Time & Storage

| Sample Type | Count | Size Each | Total | Est. Time/Sample |
|-------------|-------|-----------|-------|------------------|
| SingleAnalyte WGS | 3 | ~100 GB | 300 GB | 4-6 hours |
| SingleAnalyte WGEM | 3 | ~100 GB | 300 GB | 4-6 hours |
| TrinitySeq (merged) | 6 | ~50-80 GB | 400 GB | 3-4 hours |
| HairyTNA | 4 | ~20 GB | 80 GB | 1-2 hours |

**Total input:** ~1 TB
**Working space needed:** ~3 TB (input + intermediate + output)
**Total estimated runtime:** 50-70 hours (sequential) or 12-16 hours (parallel)

## Notes

- Use `--biscuit_chunk_size 0` to disable scatter-gather (no splitting)
- Use `--revelio_chunk_size 0` to disable Revelio scatter-gather
- Add `-resume` to continue from checkpoints if interrupted
- Monitor with `htop` and `df -h` to ensure resources are sufficient
