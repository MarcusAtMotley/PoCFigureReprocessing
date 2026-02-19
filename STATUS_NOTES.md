# Pipeline Execution Status — Feb 19, 2026 (~20:30 UTC)

## ALL PIPELINES COMPLETE

| Pipeline | Complete | Running | Waiting | Total |
|----------|----------|---------|---------|-------|
| P1 DNA SNP | **16** | 0 | 0 | **16 COMPLETE** |
| P2 DNA Meth | **13** | 0 | 0 | **13 COMPLETE** |
| P3 CNV | **16** | 0 | 0 | **16 COMPLETE** |
| P4 RNA Counts | **15** | 0 | 0 | **15 COMPLETE** |
| P5 RNA SNP | **15** | 0 | 0 | **15 COMPLETE** |

**Note:** P1 expanded to 16 samples — 3 WGEM samples now run through P1 SNP calling using markdup BAMs from P2.

## MCP Integration (added Feb 13)
- **Seqera MCP**: List/inspect Tower runs directly (replaces `tw` CLI)
- **AWS MCP (`call_aws`)**: Clean Batch job listing/description (no temp file parsing)
- **CloudWatch Log Insights MCP**: Query job logs with proper filtering/sorting (required adding `logs:StartQuery`, `logs:StopQuery`, `logs:GetQueryResults` IAM permissions)

## P1 DNA SNP — 16/16 COMPLETE

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
| HT29_01Z_1A3_1DNA (WGEM) | AWS Batch | Feb 17 ~22:54 (24.9h) |
| CoB_01W_1A3_1DNA (WGEM) | AWS Batch | Feb 18 ~16:50 (24.4h) |
| CoM_01T_1A3_1DNA (WGEM) | AWS Batch | Feb 18 ~20:45 (28.3h) |
| HT29_02N_1B3_1DNA (WGS) | Local EC2 (92GB) | Feb 19 20:29 (~26.7h) — 5,404,221 variants |

### P1 WGS OOM History (HT29_02N — 140GB sorted BAM, 2B reads)
1. **Original run**: OOM during markdup (60GB RAM, 32 threads)
2. **Resume v1**: 120GB RAM rejected — c5a.8xlarge max 64GB
3. **Resume v2**: 60GB + 8 threads — OOM during sort ("merging from 22 files and 16 in-memory blocks")
4. **Resume v3**: 60GB + 4 threads + 256M sort memory — sort succeeded but **markdup OOM'd** (hash table + 140GB page cache)
5. **v4 (GATK MarkDuplicates)**: NullPointerException — biscuit BAMs have no RG tags.
6. **v5 (piped)**: `collate|fixmate|sort|markdup`. OOM'd — hash table alone >60GB for 2B reads.
7. **v6 (2-way split)**: chr1-11/chr12-Y piped. OOM'd — chr1-11 = 63% of genome, hash ~45GB + page cache.
8. **v7 (3-way split)**: Physical files chr1-6/chr7-14/chr15-Y, delete original. OOM'd — Group A (60GB) still too big.
9. **v8 (5-way split + 2-stage)**: Separated sort and markdup into stages with page cache drops. OOM'd — all group BAMs in page cache from extraction overwhelmed cgroups v1.
10. **v9 (local EC2)**: 92GB instance, no cgroups v1 page cache limits. **SUCCEEDED** — markdup → calmd → revelio (6 chromosome groups) → bcftools mpileup+call.

### P1 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p1_dna_snp/{SAMPLE}/`
- `{SAMPLE}.sorted.bam` — alignment checkpoint (input for P2/P3)
- `{SAMPLE}.revelio.sorted.bam` — post-revelio BAM
- `{SAMPLE}.bcftools.vcf.gz` — SNP calls

## P2 DNA Methylation — 13/13 COMPLETE

### Completed P2 Samples
- 7 mTNA (BAM-mode) — Feb 11
- 3 HairyTNA (BAM-mode) — Feb 11
- HT29_01Z_1A3_1DNA (WGEM FASTQ-mode) — Feb 15 (~76.6h)
- CoB_01W_1A3_1DNA (WGEM resub BAM-mode) — Feb 17 (~7.8h)
- CoM_01T_1A3_1DNA (WGEM resub BAM-mode) — Feb 17 (~8.6h)

### P2 WGEM Failure & Resubmission History
Two of three WGEM FASTQ-mode jobs failed after completing alignment + markdup:

| Sample | Original Failure | Root Cause | Fix |
|--------|-----------------|------------|-----|
| CoB_01W | OOM (exit 137) | `biscuit pileup -@ 32` on 115GB BAM | Reduced to 16 threads |
| CoM_01T | exit code 1 | Unsorted VCF positions broke tabix | Added `bcftools sort` before tabix |

Both resubmitted as BAM-mode from sorted BAMs. Succeeded Feb 17.

### P2 Output Location
`s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p2_dna_meth/{SAMPLE}/`
- `{SAMPLE}.methylation.vcf.gz` + `.tbi`
- `{SAMPLE}.markdup.bam` + `.bai` (uploaded before pileup as checkpoint)

## P3 CNV (IchorCNA) — 16/16 COMPLETE

### Completed P3 (16)
- 9 mTNA + HairyTNA batch (Feb 12)
- 1 test job (Feb 12)
- CoM_02K_1C3_1DNA (Feb 14, 5.4h from sorted BAM)
- HT29_01Z_1A3_1DNA (Feb 16, 42 min from markdup BAM)
- CoB_02M_1C3_1DNA (Feb 17, 6.0h from sorted BAM)
- CoB_01W_1A3_1DNA (Feb 17, 42 min from markdup BAM)
- CoM_01T_1A3_1DNA (Feb 17, 48 min from markdup BAM)
- HT29_02N_1B3_1DNA (Feb 19, ~35 min from markdup BAM — auto-submitted when local markdup uploaded to S3)

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
- [x] ~~P2 CoB_01W + CoM_01T resubs~~ — SUCCEEDED Feb 17 (7.8h, 8.6h)
- [x] ~~P3 CoB_02M~~ — SUCCEEDED Feb 17 (6.0h)
- [x] ~~P3 HT29_01Z~~ — SUCCEEDED (42 min from markdup BAM)
- [x] ~~P1 CoB_01W, CoM_01T from markdup~~ — SUCCEEDED Feb 18 (24.4h, 28.3h)
- [x] ~~P3 for CoB_01W, CoM_01T~~ — SUCCEEDED Feb 17 (42 min, 48 min from markdup BAMs)
- [x] ~~Generate intergenic BED file~~ — Created `references/intergenic_regions.hg38.bed` (32,661 regions, 42.3% of genome). Uploaded to S3.
- [x] ~~P1 HT29_02N~~ — SUCCEEDED Feb 19 20:29 (local EC2, 5.4M variants). 9 OOM attempts on Batch before local success.
- [x] ~~P1 HT29_01Z~~ — SUCCEEDED Feb 17 22:54 (24.9h from markdup BAM)
- [x] ~~P3 for HT29_02N~~ — SUCCEEDED Feb 19 (~35 min from markdup BAM, auto-submitted)
- [ ] **PoC figure analysis**: All pipelines complete. P1/P2/P3/P4/P5 data ready on S3.

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
- **samtools markdup used for all samples** (not GATK/Picard) — biscuit BAMs lack RG tags, and samtools markdup produces equivalent results.

## Lessons Learned
- **Always test with a single sample before batch submission.** Multiple P2 and P3 failures wasted compute.
- The biscuit index directory's FASTA (`biscuit_reference_genome/`) is a Fusion symlink placeholder (90 bytes). Use regular reference.
- `biscuit pileup` does NOT accept `-q` flag. Base quality = `-b`, mapping quality = `-m`.
- AWS Batch reuses EC2 instances — `/scratch` can have stale files. Always clean at script start.
- cnvpytor has too many dependency issues (pkg_resources, numpy 2.0, pickle). GATK fails on sparse BAMs. **IchorCNA is the right tool for mixed-coverage CNV**.
- IchorCNA conda package (`r-ichorcna`) does NOT include `runIchorCNA.R` CLI script. Call `run_ichorCNA()` R function directly. The readCounter binary is in `hmmcopy` (not `hmmcopy-utils`).
- Set `genomeBuild='hg38'`, `genomeStyle='UCSC'` for chr-prefixed BAMs in IchorCNA.
- **WGS 140GB BAM markdup OOM saga (9 attempts)**: Fundamentally impossible on 60GB container — hash table for 2B reads is ~72GB, and cgroups v1 counts file page cache against container memory limit. Splitting by chromosome doesn't help enough because page cache from extracted group BAMs still overwhelms. **Fix: run locally on EC2 instance with 92GB RAM** (no cgroups v1 overhead).
- **GATK/Picard MarkDuplicates requires RG tags** — biscuit align doesn't add them. Use samtools markdup instead.
- **Biscuit pileup OOMs at 32 threads** on large WGEM BAMs. Fix: use 16 threads.
- **Biscuit pileup can produce unsorted VCF positions.** Fix: `bcftools sort` before tabix.
- **Upload checkpoints before expensive steps.** Markdup BAM should be uploaded before pileup, sorted BAM before markdup.
- **Biscuit auto-detects gzip** by magic bytes, not file extension.
- **WGEM alignment is ~2x slower per read than WGS** through biscuit — bisulfite 3-letter alignment.
- **MCP tools**: Seqera, AWS `call_aws`, and CloudWatch Log Insights MCPs greatly simplify monitoring.
- **cgroups v1 page cache**: Container memory limits include file page cache. Large intermediate BAMs (140GB+) consume page cache that competes with application memory, causing OOM even when the app alone fits in memory. Solution: pipe operations to avoid intermediate files.
