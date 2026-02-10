# Full Pipeline Execution on VM

Run P1 (DNA SNP), P2 (DNA Methylation), P3 (CNV), and P5 (RNA SNP) on a large VM.
P4 (RNA Counts) is already complete.

## Pipeline Overview

| Pipeline | Input | Output | Status |
|----------|-------|--------|--------|
| P1: DNA SNP | DNA FASTQs | VCF | To run |
| P2: DNA Methylation | DNA FASTQs | VCF | To run |
| P3: CNV | DNA FASTQs | CSV | To run |
| P4: RNA Counts | RNA FASTQs | Counts | ✅ DONE |
| P5: RNA SNP | RNA BAMs (from P4) | VCF | To run |

## VM Requirements

- **Instance:** i3.8xlarge (32 vCPU, 244 GB RAM, 7.6 TB NVMe)
- **Region:** us-east-2
- **Storage needed:** ~4 TB (inputs + working + outputs)

## Quick Start

```bash
# On VM
git clone https://github.com/MarcusAtMotley/PoCFigureReprocessing.git
cd PoCFigureReprocessing

./scripts/vm_setup.sh
# Log out/in for docker group, then:
aws configure

./scripts/vm_download_all.sh      # Download everything
./scripts/vm_run_all_pipelines.sh # Run P1, P2, P3, P5
./scripts/vm_upload_results.sh    # Upload to S3
```

## Sample Inventory

### DNA Samples (P1, P2, P3)

| Sample | Type | Pipeline | Cell Line |
|--------|------|----------|-----------|
| CoB_02M_1C3_1DNA | SingleAnalyte WGS | P1 | CoB |
| CoM_02K_1C3_1DNA | SingleAnalyte WGS | P1 | CoM |
| HT29_02N_1B3_1DNA | SingleAnalyte WGS | P1 | HT29 |
| CoB_01W_1A3_1DNA | SingleAnalyte WGEM | P2 | CoB |
| CoM_01T_1A3_1DNA | SingleAnalyte WGEM | P2 | CoM |
| HT29_01Z_1A3_1DNA | SingleAnalyte WGEM | P2 | HT29 |
| CoB_08L_3A2_DNA-EM | TrinitySeq | P1\|P2 | CoB |
| CoM_08M_3A2_DNA-EM | TrinitySeq | P1\|P2 | CoM |
| CoB_08R_3A2_TNA-mRT-EM | TrinitySeq mTNA | P1\|P2 | CoB |
| CoM_08S_3A2_TNA-mRT-EM | TrinitySeq mTNA | P1\|P2 | CoM |
| CoB_08X_3A2_TNA-RT-EM | TrinitySeq TNA-RT | P1\|P2 | CoB |
| CoM_08Y_3A2_TNA-RT-EM | TrinitySeq TNA-RT | P1\|P2 | CoM |
| HT29_21S_3A2_bsTNA-HP | HairyTNA | P1\|P2 | HT29 |
| HT29_21T_3A2_bsTNA-HP | HairyTNA | P1\|P2 | HT29 |
| HT29_21W_3A2_bsDNA-HP | HairyTNA | P1\|P2 | HT29 |
| HT29_21X_3A2_bsDNA-HP | HairyTNA | P1\|P2 | HT29 |

### RNA Samples (P5) - Using BAMs from completed P4

| Sample | Type | Cell Line |
|--------|------|-----------|
| CoB_02W_1A3_1RNA | SingleAnalyte | CoB |
| CoM_02V_1A3_1RNA | SingleAnalyte | CoM |
| HT29_02T_1A3_1RNA | SingleAnalyte | HT29 |
| CoB_08R_3A2_TNA-mRT-EM | mTNA | CoB |
| CoM_08S_3A2_TNA-mRT-EM | mTNA | CoM |
| HT29_21S_3A2_bsTNA-HP | HairyTNA | HT29 |
| HT29_21T_3A2_bsTNA-HP | HairyTNA | HT29 |
| HT29_21U_3A2_bsRNA-HP | HairyTNA | HT29 |
| HT29_21V_3A2_bsRNA-HP | HairyTNA | HT29 |

## Directory Structure on VM

```
/data/
├── references/
│   ├── GRCh38_full_analysis_set_plus_decoy_hla.fa
│   ├── GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
│   └── biscuit_index/
├── fastq/                    # DNA FASTQs for P1/P2/P3
├── rna_bams/                 # RNA BAMs from P4 for P5
├── results/
│   ├── p1_dna_snp/
│   ├── p2_dna_meth/
│   ├── p3_cnv/
│   └── p5_rna_snp/
└── PoCFigureReprocessing/    # Pipeline code
```

## Estimated Runtime

| Pipeline | Samples | Est. Time |
|----------|---------|-----------|
| P1 DNA SNP | 16 | 24-36 hours |
| P2 DNA Meth | 16 | 12-18 hours |
| P3 CNV | 16 | 8-12 hours |
| P5 RNA SNP | 9 | 12-18 hours |
| **Total** | - | **~60-80 hours** |

Running in parallel where possible can reduce this significantly.

## Notes

- Use `-resume` to continue from checkpoints if interrupted
- Monitor disk space with `df -h /data`
- Monitor progress with `htop` and `tail -f .nextflow.log`
- Results upload to: `s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/`
