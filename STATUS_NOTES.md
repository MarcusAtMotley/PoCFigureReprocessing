# Pipeline Execution Status — Feb 11, 2026

## P1 DNA SNP — Nearly Complete

**10/13 samples done.** 3 WGS SingleAnalyte still aligning on AWS Batch (~31+ hours, full 30x WGS through biscuit).

### Completed P1 Samples
| Sample | How | sorted.bam on S3? |
|--------|-----|-------------------|
| CoB_08L_3A2_DNA-EM | Local EC2 | Yes (uploaded manually) |
| CoB_08R_3A2_TNA-mRT-EM | Local EC2 | Yes (uploaded manually) |
| CoM_08M_3A2_DNA-EM | Local EC2 | Yes (uploaded manually) |
| CoB_08X_3A2_TNA-RT-EM | AWS Batch | Yes |
| CoM_08S_3A2_TNA-mRT-EM | AWS Batch | Yes |
| CoM_08Y_3A2_TNA-RT-EM | AWS Batch | Yes |
| HT29_21S_3A2_bsTNA-HP | AWS Batch | Yes |
| HT29_21T_3A2_bsTNA-HP | AWS Batch | Yes |
| HT29_21W_3A2_bsDNA-HP | AWS Batch | Yes |
| HT29_21X_3A2_bsDNA-HP | AWS Batch | Yes |

### Still Running P1 (AWS Batch, queue: C4_QUEUE)
| Sample | Job ID | Stage |
|--------|--------|-------|
| CoB_02M_1C3_1DNA (WGS) | c0bf1e59-77cb-42c9-893c-4b7b45b50660 | Biscuit alignment |
| CoM_02K_1C3_1DNA (WGS) | d32c6149-753c-457c-b462-c07440885ba3 | Biscuit alignment |
| HT29_02N_1B3_1DNA (WGS) | 752c291d-1499-4a34-afbf-28ae749a8457 | Biscuit alignment |

### P1 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p1_dna_snp/{SAMPLE}/`
- `{SAMPLE}.sorted.bam` — alignment checkpoint (input for P2/P3)
- `{SAMPLE}.revelio.sorted.bam` — post-revelio BAM
- `{SAMPLE}.bcftools.vcf.gz` — SNP calls

## P2 DNA Methylation — SUBMITTED (Feb 11, ~18:30 UTC)

**All 13 jobs submitted to AWS Batch.** Job queue: C4_QUEUE.
- 10 BAM-mode samples reuse P1 sorted BAMs (all on S3)
- 3 WGEM FASTQ-mode samples align from scratch (no P1 dependency)

### How to Submit P2
```bash
# Docker image must be rebuilt first (adds cnvpytor, P2/P3 scripts)
docker build -f Dockerfile.batch -t poc-dna-pipeline:latest .
docker tag poc-dna-pipeline:latest 891377403536.dkr.ecr.us-east-2.amazonaws.com/poc-dna-pipeline:latest
docker push 891377403536.dkr.ecr.us-east-2.amazonaws.com/poc-dna-pipeline:latest

# Submit P2 jobs
./scripts/submit_p2_p3_jobs.sh p2

# Or dry run first
./scripts/submit_p2_p3_jobs.sh p2 --dry-run
```

### P2 Samplesheet: `samplesheets/p2_batch.csv`
- 10 BAM-mode (sorted.bam → markdup → biscuit pileup)
- 3 FASTQ-mode WGEM (align → markdup → biscuit pileup)

### P2 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p2_dna_meth/{SAMPLE}/`
- `{SAMPLE}.methylation.vcf.gz` + `.tbi`
- `{SAMPLE}.markdup.bam` + `.bai` (saved for P3 reuse)

## P3 CNV — After P2

**16 samples total.** Depends on P2 markdup BAMs for 13 samples; 3 WGS use P1 sorted BAMs directly.

### How to Submit P3
```bash
# After P2 completes:
./scripts/submit_p2_p3_jobs.sh p3
```

### P3 Samplesheet: `samplesheets/p3_batch.csv`
- 13 markdup BAMs from P2 output
- 3 WGS sorted BAMs from P1 (will run markdup internally)

### P3 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p3_cnv/{SAMPLE}/`

## Key Infrastructure
- **Job queue:** `C4_QUEUE`
- **Job definition:** `poc-dna-pipeline` (ECR: `891377403536.dkr.ecr.us-east-2.amazonaws.com/poc-dna-pipeline:latest`)
- **Docker image:** `Dockerfile.batch` — includes biscuit, samtools, bcftools, cnvpytor, bedtools, revelio

## Decision Log
- Running all samples at **full depth** (no downsampling). Can downsample at analysis stage later if needed for fair comparisons.
- P2 uses **markdup BAM** (pre-revelio) — revelio would destroy methylation signal.
- P3 also uses **markdup BAM** — read depth analysis needs clean BAM.
- P3 CNV calling is **raw** (no intergenic filtering). Intergenic masking for RNA bleedthrough will be done post-hoc — it's just a `bedtools intersect` and takes seconds. Need to generate the intergenic BED file from GTF first.
- cnvpytor has a `pkg_resources` import error in the Docker image — needs fix before P3 submission (add `pip install setuptools` to Dockerfile).
