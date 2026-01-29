#!/bin/bash
export TOWER_ACCESS_TOKEN='eyJ0aWQiOiAxMzIxM30uNTZjMDEyZmI3ZjkwMjE0Yjk4ZWE3NGFlNDkyM2IwOGJiYTk0M2JmMg=='

cd /Users/marcus/PycharmProjects/MotleyDataFigures/poc_reprocessing

nextflow run main.nf \
  -profile awsbatch \
  --input "s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/samplesheets/p4_extended_samples.csv" \
  --outdir "s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS" \
  --genome_fasta "s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --genome_fai "s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" \
  --annotation_gtf "s3://motleybio/Resources/GTF/Homo_sapiens.GRCh38.112.chr_label.gtf" \
  --star_index "s3://motleybio/Resources/Star_Genome/" \
  --aws_queue "TowerForge-3kv8myJytVMvWpM0JA6lcn-work" \
  -w s3://motleybio-nf-work/scratch/local_p4_extended \
  -resume
