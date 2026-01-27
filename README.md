# PoC Data Reprocessing

Data reprocessing plan and manifests for generating comparable PoC figures across Single Analyte, mTNA, and HairyTNA library types.

## Quick Start

1. Read `DATA_REPROCESSING_PLAN.md` for full context
2. Use `sample_manifest.csv` for raw FASTQ inventory
3. Use `pipeline_ready_manifest.csv` for pipeline-specific inputs

## Files

| File | Description |
|------|-------------|
| `DATA_REPROCESSING_PLAN.md` | Comprehensive plan with data sources, pipelines, terminology |
| `sample_manifest.csv` | Raw FASTQ inventory with S3 paths |
| `pipeline_ready_manifest.csv` | Pipeline-ready FASTQs (post-split intermediates) |
| `processing_pipelines.md` | Pipeline documentation |

## Cell Lines

- **CoB** - Single analyte + mTNA
- **CoM** - Single analyte + mTNA
- **HT29** - Single analyte + HairyTNA (bisulfite)

## Data Locations

- Single Analyte: `s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/`
- mTNA (mot26): `s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley26/results/`
- HairyTNA (mot34): `s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/`

## Pipelines

1. **P1: DNA SNP** - Biscuit → Revelio → LoFreq
2. **P2: DNA Methylation** - Biscuit → Biscuit Pileup
3. **P3: CNV** - Hisat2 → CNAnator → Intergenic Filter
4. **P4: RNA Counts** - STAR → FeatureCounts
5. **P5: RNA SNP** - STAR → Revelio → LoFreq
