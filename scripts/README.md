# Pipeline Scripts

## Current Scripts

### Main Pipelines

| Script | Description |
|--------|-------------|
| `simple_p5_rna_snp.sh` | P5 RNA SNP pipeline: BAM → calmd → Revelio → BCFtools → VCF |
| `simple_dna_pipeline.sh` | P1/P2/P3 DNA pipeline: FASTQ → Biscuit/BWA → MarkDup → Revelio → BCFtools/Biscuit/CNVnator |

### Runners (use these to execute)

| Script | Description |
|--------|-------------|
| `run_p5_trinityseq.sh` | Run P5 on 12 TrinitySeq RNA samples |
| `run_p5_singleanalyte.sh` | Run P5 on 3 SingleAnalyte RNA samples (downsampled to 30M) |
| `run_p5_all.sh` | Run P5 on all 15 RNA samples |

### Utilities

| Script | Description |
|--------|-------------|
| `downsample_bam.sh` | Downsample BAM to target read count |

## Quick Start

```bash
# P5 RNA SNP - start with TrinitySeq (smaller, faster)
./run_p5_trinityseq.sh

# Then SingleAnalyte (larger, downsampled)
./run_p5_singleanalyte.sh

# Or all at once
./run_p5_all.sh
```

## Key Design Decisions

1. **BCFtools** replaces LoFreq (50-100x faster)
2. **Revelio applied to ALL samples** for pipeline consistency
3. **SingleAnalyte downsampled** to match TrinitySeq depths (RNA: 30M, DNA: 140M)
4. **Quality thresholds**: RNA q13/Q13, DNA q20/Q20

## Archived Scripts

Legacy scripts moved to `archive/` - kept for reference but not actively used.
