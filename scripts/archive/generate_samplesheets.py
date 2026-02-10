#!/usr/bin/env python3
"""
Generate pipeline-specific samplesheets from the pipeline_ready_manifest.csv.

This script transforms the manifest into per-pipeline CSV files in nf-core format:
    sample,fastq_1,fastq_2,single_end,pipeline,cell_line,assay_category

Handles lane merging: L001/L002 samples are combined into single entries with
comma-separated FASTQ paths that will be merged by the pipeline.

Usage:
    python generate_samplesheets.py \
        --manifest ../pipeline_ready_manifest.csv \
        --output-dir ../samplesheets

Author: Motley Bio
"""

import argparse
import csv
import re
from pathlib import Path
from collections import defaultdict


def parse_manifest(manifest_path: Path) -> list[dict]:
    """Parse the pipeline_ready_manifest.csv file."""
    samples = []

    with open(manifest_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Skip comment lines
            if row.get('sample_id', '').startswith('#'):
                continue
            # Skip empty rows
            if not row.get('sample_id'):
                continue

            samples.append({
                'sample_id': row['sample_id'],
                'cell_line': row['cell_line'],
                'assay_category': row['assay_category'],
                'target_channel': row['target_channel'],
                'pipeline': row['pipeline'],
                'read_type': row['read_type'],
                'fastq_path': row['fastq_path'],
                'fastq_r2_path': row.get('fastq_r2_path', ''),
                'needs_processing': row.get('needs_processing', 'ready')
            })

    return samples


def get_base_sample_id(sample_id: str) -> str:
    """Remove lane suffix (L001, L002, etc.) from sample ID."""
    # Match patterns like _L001, _L002 at the end
    return re.sub(r'_L00\d$', '', sample_id)


def merge_lanes(samples: list[dict]) -> list[dict]:
    """Merge samples that differ only by lane number."""

    # Group samples by (base_sample_id, pipeline)
    grouped = defaultdict(list)
    for sample in samples:
        base_id = get_base_sample_id(sample['sample_id'])
        key = (base_id, sample['pipeline'])
        grouped[key].append(sample)

    merged_samples = []
    for (base_id, pipeline), group in grouped.items():
        if len(group) == 1:
            # No merging needed, use original sample_id
            merged_samples.append(group[0])
        else:
            # Multiple lanes - merge FASTQs
            # Sort by lane to ensure consistent ordering
            group.sort(key=lambda x: x['sample_id'])

            # Combine FASTQ paths with semicolons (Nextflow-friendly)
            r1_paths = ';'.join(s['fastq_path'] for s in group)
            r2_paths = ';'.join(s['fastq_r2_path'] for s in group if s['fastq_r2_path'])

            # Use first sample as template, update with merged info
            merged = group[0].copy()
            merged['sample_id'] = base_id  # Remove lane suffix
            merged['fastq_path'] = r1_paths
            merged['fastq_r2_path'] = r2_paths if r2_paths else ''
            merged['lanes_merged'] = len(group)

            merged_samples.append(merged)

    return merged_samples


def generate_samplesheet(samples: list[dict], pipeline: str, output_path: Path):
    """Generate a samplesheet for a specific pipeline."""

    # Filter samples for this pipeline
    pipeline_samples = [s for s in samples if s['pipeline'] == pipeline]

    if not pipeline_samples:
        print(f"  No samples for {pipeline}, skipping")
        return 0

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'fastq_1', 'fastq_2', 'single_end', 'pipeline', 'cell_line', 'assay_category'])

        for sample in pipeline_samples:
            single_end = sample['read_type'] == 'SE'
            fastq_2 = '' if single_end else sample['fastq_r2_path']

            writer.writerow([
                sample['sample_id'],
                sample['fastq_path'],
                fastq_2,
                str(single_end).lower(),
                sample['pipeline'],
                sample['cell_line'],
                sample['assay_category']
            ])

    print(f"  {output_path.name}: {len(pipeline_samples)} samples")
    return len(pipeline_samples)


def generate_combined_samplesheet(samples: list[dict], output_path: Path):
    """Generate a combined samplesheet with all samples."""

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'fastq_1', 'fastq_2', 'single_end', 'pipeline', 'cell_line', 'assay_category'])

        for sample in samples:
            single_end = sample['read_type'] == 'SE'
            fastq_2 = '' if single_end else sample['fastq_r2_path']

            writer.writerow([
                sample['sample_id'],
                sample['fastq_path'],
                fastq_2,
                str(single_end).lower(),
                sample['pipeline'],
                sample['cell_line'],
                sample['assay_category']
            ])

    print(f"  {output_path.name}: {len(samples)} samples (all pipelines)")


def generate_p0_samplesheet(samples: list[dict], output_path: Path):
    """Generate P0 (TrimGalore) samplesheet for samples needing preprocessing."""

    # Get unique samples that need TrimGalore (single analyte samples)
    p0_samples = []
    seen_samples = set()

    for sample in samples:
        base_id = get_base_sample_id(sample['sample_id'])
        if sample['needs_processing'] == 'TrimGalore' and base_id not in seen_samples:
            seen_samples.add(base_id)
            p0_samples.append(sample)

    if not p0_samples:
        print(f"  No samples need TrimGalore preprocessing")
        return

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'fastq_1', 'fastq_2', 'single_end', 'pipeline', 'cell_line', 'assay_category'])

        for sample in p0_samples:
            single_end = sample['read_type'] == 'SE'
            fastq_2 = '' if single_end else sample['fastq_r2_path']

            writer.writerow([
                sample['sample_id'],
                sample['fastq_path'],
                fastq_2,
                str(single_end).lower(),
                'P0_TrimGalore',
                sample['cell_line'],
                sample['assay_category']
            ])

    print(f"  {output_path.name}: {len(p0_samples)} unique samples")


def main():
    parser = argparse.ArgumentParser(
        description='Generate pipeline-specific samplesheets from manifest'
    )
    parser.add_argument(
        '--manifest', '-m',
        type=Path,
        default=Path('../pipeline_ready_manifest.csv'),
        help='Path to pipeline_ready_manifest.csv'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=Path('../samplesheets'),
        help='Output directory for samplesheets'
    )
    parser.add_argument(
        '--no-merge-lanes',
        action='store_true',
        help='Do not merge L001/L002 lanes (keep separate)'
    )

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Parse manifest
    print(f"Reading manifest: {args.manifest}")
    samples = parse_manifest(args.manifest)
    print(f"Found {len(samples)} sample-pipeline rows in manifest")

    # Merge lanes unless disabled
    if not args.no_merge_lanes:
        samples = merge_lanes(samples)
        print(f"After lane merging: {len(samples)} sample-pipeline combinations\n")
    else:
        print(f"Lane merging disabled: {len(samples)} sample-pipeline combinations\n")

    # Generate samplesheets
    print("Generating samplesheets:")

    # P0: TrimGalore (unique samples needing preprocessing)
    generate_p0_samplesheet(samples, args.output_dir / 'p0_trimgalore.csv')

    # Per-pipeline samplesheets
    total_runs = 0
    pipelines = ['P1_DNA_SNP', 'P2_DNA_Meth', 'P3_CNV', 'P4_RNA_Counts', 'P5_RNA_SNP']
    for pipeline in pipelines:
        output_file = args.output_dir / f"{pipeline.lower()}.csv"
        count = generate_samplesheet(samples, pipeline, output_file)
        total_runs += count

    # Combined samplesheet (all samples, ready only)
    ready_samples = [s for s in samples if s['needs_processing'] == 'ready']
    generate_combined_samplesheet(ready_samples, args.output_dir / 'all_ready_samples.csv')

    print(f"\nTotal pipeline runs: {total_runs}")
    print("Done!")


if __name__ == '__main__':
    main()
