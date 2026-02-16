# Pipeline Execution Status — Feb 16, 2026 (~20:00 UTC)

## Overall Progress

| Pipeline | Complete | Running | Waiting | Total |
|----------|----------|---------|---------|-------|
| P1 DNA SNP | 12 | 1 (HT29_02N resume v3) | 0 | 13 |
| P2 DNA Meth | 11 | 2 (CoB_01W + CoM_01T resubs from sorted BAM) | 0 | 13 |
| P3 CNV | 11 | 1 (CoB_02M) | 4 | 16 |
| P4 RNA Counts | **15** | 0 | 0 | **15 COMPLETE** |
| P5 RNA SNP | **15** | 0 | 0 | **15 COMPLETE** |

## MCP Integration (added Feb 13)
- **Seqera MCP**: List/inspect Tower runs directly (replaces `tw` CLI)
- **AWS MCP (`call_aws`)**: Clean Batch job listing/description (no temp file parsing)
- **CloudWatch Log Insights MCP**: Query job logs with proper filtering/sorting (required adding `logs:StartQuery`, `logs:StopQuery`, `logs:GetQueryResults` IAM permissions)

## P1 DNA SNP — 12/13 Done

### Completed P1 Samples
| Sample | How | Finished |
|--------|-----|----------|
| CoB_08L_3A2_DNA-EM | Local EC2 | ~Feb 8 |
| CoB_08R_3A2_TNA-mRT-EM | Local EC2 | ~Feb 8 |
| CoM_08M_3A2_DNA-EM | Local EC2 | ~Feb 8 |
| CoB_08X_3A2_TNA-RT-EM | AWS Batch | Feb 10 |
| CoM_08S_3A2_TNA-mRT-EM | AWS Batch | Feb 11 |
| CoM_08Y_3A2_TNA-RT-EM | AWS Batch | Feb 11 |
| HT29_21S_3A2_bsTNA-HP | AWS Batch | Feb 11 |
| HT29_21T_3A2_bsTNA-HP | AWS Batch | Feb 11 |
| HT29_21W_3A2_bsDNA-HP | AWS Batch | Feb 11 |
| HT29_21X_3A2_bsDNA-HP | AWS Batch | Feb 11 |
| CoM_02K_1C3_1DNA (WGS) | AWS Batch | Feb 13 19:36 (84.3h) |
| CoB_02M_1C3_1DNA (WGS) | AWS Batch | Feb 14 ~01:00 (90.0h) |

### Still Running P1
| Sample | Job ID | Notes |
|--------|--------|-------|
| HT29_02N_1B3_1DNA (WGS) | ee3ff705 (v3) | Resume from sorted BAM. **4 threads + 256M sort memory.** Previous attempts: OOM at 32 threads, OOM at 8 threads. |

### P1 WGS OOM History (HT29_02N — 150GB sorted BAM)
1. **Original run**: OOM during markdup (60GB RAM, 32 threads)
2. **Resume v1**: 120GB RAM rejected — c5a.8xlarge max 64GB
3. **Resume v2**: 60GB + 8 threads — OOM during sort ("merging from 22 files and 16 in-memory blocks")
4. **Resume v3 (current)**: 60GB + 4 threads + explicit `-m 256M` sort memory. Also deletes input BAM after collate to free page cache.

### P1 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p1_dna_snp/{SAMPLE}/`
- `{SAMPLE}.sorted.bam` — alignment checkpoint (input for P2/P3)
- `{SAMPLE}.revelio.sorted.bam` — post-revelio BAM
- `{SAMPLE}.bcftools.vcf.gz` — SNP calls

## P2 DNA Methylation — 11/13 Complete

### Completed P2 Samples
- 7 mTNA (BAM-mode) — Feb 11
- 3 HairyTNA (BAM-mode) — Feb 11
- HT29_01Z_1A3_1DNA (WGEM FASTQ-mode) — Feb 15 (~76.6h)

### P2 WGEM Failure & Resubmission
Two of three WGEM FASTQ-mode jobs failed after completing alignment + markdup:

| Sample | Original Failure | Root Cause | Fix |
|--------|-----------------|------------|-----|
| CoB_01W | OOM (exit 137) | `biscuit pileup -@ 32` on 115GB BAM | Reduced to 16 threads |
| CoM_01T | exit code 1 | Unsorted VCF positions broke tabix | Added `bcftools sort` before tabix |

**Both have sorted BAMs on S3** from alignment checkpoint. Resubmitted as BAM-mode.

### P2 Script Fixes (Feb 16)
1. Markdup: 4 threads + `-m 256M` sort memory (was 32 threads, default 768M)
2. Biscuit pileup: 16 threads (was 32 — OOM'd CoB_01W)
3. Added `bcftools sort` before tabix (fixes unsorted VCF from CoM_01T)
4. Markdup BAM uploaded as checkpoint BEFORE pileup (saves markdup work if pileup fails)

### Still Running P2
| Sample | Job ID | Mode | Notes |
|--------|--------|------|-------|
| CoB_01W_1A3_1DNA | ca13868f | BAM (resub) | From sorted BAM. 4-thread markdup, 16-thread pileup. |
| CoM_01T_1A3_1DNA | 37d6658a | BAM (resub) | From sorted BAM. Same fixes + VCF sort. |

### P2 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p2_dna_meth/{SAMPLE}/`
- `{SAMPLE}.methylation.vcf.gz` + `.tbi`
- `{SAMPLE}.markdup.bam` + `.bai` (now uploaded before pileup as checkpoint)

## P3 CNV (IchorCNA) — 11/16 Complete, 1 Running

### Completed P3 (11)
- 9 mTNA + HairyTNA batch (Feb 12)
- 1 test job (Feb 12)
- CoM_02K_1C3_1DNA (Feb 14, 5.4h from sorted BAM)

### Now Running P3
| Sample | Job ID | Notes |
|--------|--------|-------|
| CoB_02M_1C3_1DNA | 766940d3 | From sorted BAM → markdup → IchorCNA. RUNNABLE. |

### Waiting on Dependencies (4 samples)
| Sample | Waiting For |
|--------|-------------|
| CoB_01W_1A3_1DNA | P2 resub completion (markdup BAM) |
| CoM_01T_1A3_1DNA | P2 resub completion (markdup BAM) |
| HT29_01Z_1A3_1DNA | Ready — P2 succeeded, markdup BAM on S3 |
| HT29_02N_1B3_1DNA | P1 resume v3 completion (sorted BAM → markdup) |

**Note:** HT29_01Z P3 can be submitted now! Markdup BAM is at `s3://motleybio/.../p2_dna_meth/HT29_01Z_1A3_1DNA/HT29_01Z_1A3_1DNA.markdup.bam`

### P3 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p3_cnv/{SAMPLE}/`

## P4 RNA Counts — 15/15 COMPLETE

All 15 samples done. 9 original (Seqera runs) + 6 extended (BAM-input run 4m7zbfGmDqIFhR).

## P5 RNA SNP — 15/15 COMPLETE

Seqera run 4m7zbfGmDqIFhR SUCCEEDED (completed Feb 13 14:00 UTC). All 9 original + 6 extended LoFreq done.

### Seqera Run: P4_P5_Extended_BAM_input
- **Run ID:** 4m7zbfGmDqIFhR
- **Status:** SUCCEEDED — completed Feb 13 14:00 UTC

## TODO
- [x] ~~P4/P5 extended RNA samples~~ — 15/15 complete
- [x] ~~P3 for CoM_02K~~ — SUCCEEDED (5.4h)
- [x] ~~P2 script fixes~~ — Reduced threads, VCF sort, checkpoint upload
- [ ] **P1 HT29_02N resume v3**: Running (job ee3ff705). 4 threads + 256M sort memory.
- [ ] **P2 CoB_01W + CoM_01T resubs**: Running (jobs ca13868f, 37d6658a). From sorted BAMs.
- [ ] **P3 CoB_02M**: Running (job 766940d3). From sorted BAM.
- [ ] **P3 HT29_01Z**: Ready to submit now (markdup BAM on S3 from P2 success)
- [ ] **P3 for CoB_01W, CoM_01T**: Submit when P2 resubs complete
- [ ] **P3 for HT29_02N**: Submit when P1 resume v3 completes
- [ ] **P1 for WGEM**: Submit 3 WGEM samples through P1 SNP calling using markdup BAMs from P2 (BAM_TYPE=markdup). Samples: CoB_01W, CoM_01T, HT29_01Z.
- [ ] Generate intergenic BED file from GTF for post-hoc RNA bleedthrough masking
- [ ] Commit all changes to git

## Analysis Tools (for PoC figure generation)
When pipelines complete, these Claude Code skills are available for analysis:
- **pysam** — BAM/VCF inspection and stats extraction
- **deeptools** — NGS QC (BAM correlation, PCA, fingerprinting across assay types)
- **pydeseq2** — Differential expression for RNA counts comparison (P4)
- **matplotlib / seaborn / plotly** — Visualization for PoC figures
- **scientific-visualization** — Publication-ready figure formatting
- **statistical-analysis** — Assay-vs-assay statistical comparisons
- **polars** — Fast DataFrame processing for count matrices and variant tables
- **cosmic-database** — Validate HT29 SNP calls against known somatic mutations
- **ensembl-database** — Gene/transcript annotation lookups
- **xlsx / pdf / pptx** — Output deliverables

## Key Infrastructure
- **Job queue:** `C4_QUEUE` (c5a.8xlarge: 32 vCPU, 64GB max)
- **Job definition:** `poc-dna-pipeline` (ECR: `891377403536.dkr.ecr.us-east-2.amazonaws.com/poc-dna-pipeline:latest`)
- **Docker image:** `Dockerfile.batch` — biscuit, samtools, bcftools, gatk4, bedtools, r-ichorcna, hmmcopy, revelio
- **Slack status script:** `scripts/slack_status_check.sh`

## Decision Log
- Running all samples at **full depth** (no downsampling). Can downsample at analysis stage later.
- P2 uses **markdup BAM** (pre-revelio) — revelio would destroy methylation signal.
- P3 uses **markdup BAM** via IchorCNA — works across all coverage levels consistently.
- P3 uses **IchorCNA for all samples** (not GATK) — consistent caller for PoC comparison across assay types, handles ultra-low-pass HairyTNA.
- P3 CNV calling is **raw** (no intergenic filtering). Intergenic masking for RNA bleedthrough will be done post-hoc with `bedtools intersect`.
- IchorCNA outputs are segment-based (BED-compatible) — downstream filtering works same as GATK output.

## Lessons Learned
- **Always test with a single sample before batch submission.** Multiple P2 and P3 failures wasted compute.
- The biscuit index directory's FASTA (`biscuit_reference_genome/`) is a Fusion symlink placeholder (90 bytes). Use regular reference.
- `biscuit pileup` does NOT accept `-q` flag. Base quality = `-b`, mapping quality = `-m`.
- AWS Batch reuses EC2 instances — `/scratch` can have stale files. Always clean at script start.
- cnvpytor has too many dependency issues (pkg_resources, numpy 2.0, pickle). GATK fails on sparse BAMs. **IchorCNA is the right tool for mixed-coverage CNV**.
- IchorCNA conda package (`r-ichorcna`) does NOT include `runIchorCNA.R` CLI script. Call `run_ichorCNA()` R function directly. The readCounter binary is in `hmmcopy` (not `hmmcopy-utils`).
- Set `genomeBuild='hg38'`, `genomeStyle='UCSC'` for chr-prefixed BAMs in IchorCNA.
- **WGS BAMs >140GB OOM during markdup with 60GB RAM.** Even 8 sort threads can OOM on 150GB BAMs. Fix: 4 threads + explicit `-m 256M` per thread. C4_QUEUE uses c5a.8xlarge (64GB max).
- **Biscuit pileup OOMs at 32 threads** on large WGEM BAMs. Fix: use 16 threads.
- **Biscuit pileup can produce unsorted VCF positions.** Fix: `bcftools sort` before tabix.
- **Upload checkpoints before expensive steps.** Markdup BAM should be uploaded before pileup, sorted BAM before markdup.
- **Biscuit auto-detects gzip** by magic bytes, not file extension.
- **WGEM alignment is ~2x slower per read than WGS** through biscuit — bisulfite 3-letter alignment.
- **MCP tools**: Seqera, AWS `call_aws`, and CloudWatch Log Insights MCPs greatly simplify monitoring.
