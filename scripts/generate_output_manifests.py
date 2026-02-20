#!/usr/bin/env python3
"""
Generate output manifest CSVs for all 5 pipelines (P1-P5).

Reads sample metadata from sample_manifest.csv, checks S3 for actual output files,
and writes per-pipeline manifest CSVs with full S3 paths and completion status.
"""

import csv
import subprocess
import os
import sys
from collections import defaultdict

# ── Configuration ──────────────────────────────────────────────────────────────

S3_BASE = "s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
REGION = "us-east-2"
MANIFEST_PATH = "/home/ubuntu/PoCFigureReprocessing/sample_manifest.csv"
OUTPUT_DIR = "/home/ubuntu/PoCFigureReprocessing/output_manifests"

# Sample lists per pipeline (as specified by user)
P1_SAMPLES = [
    "CoB_08L_3A2_DNA-EM", "CoB_08R_3A2_TNA-mRT-EM", "CoB_08X_3A2_TNA-RT-EM",
    "CoM_08M_3A2_DNA-EM", "CoM_08S_3A2_TNA-mRT-EM", "CoM_08Y_3A2_TNA-RT-EM",
    "HT29_21S_3A2_bsTNA-HP", "HT29_21T_3A2_bsTNA-HP",
    "HT29_21W_3A2_bsDNA-HP", "HT29_21X_3A2_bsDNA-HP",
    "CoB_02M_1C3_1DNA", "CoM_02K_1C3_1DNA", "HT29_02N_1B3_1DNA",
    "CoB_01W_1A3_1DNA", "CoM_01T_1A3_1DNA", "HT29_01Z_1A3_1DNA",
]

P2_SAMPLES = [
    "CoB_08L_3A2_DNA-EM", "CoB_08R_3A2_TNA-mRT-EM", "CoB_08X_3A2_TNA-RT-EM",
    "CoM_08M_3A2_DNA-EM", "CoM_08S_3A2_TNA-mRT-EM", "CoM_08Y_3A2_TNA-RT-EM",
    "HT29_21S_3A2_bsTNA-HP", "HT29_21T_3A2_bsTNA-HP",
    "HT29_21W_3A2_bsDNA-HP", "HT29_21X_3A2_bsDNA-HP",
    "CoB_01W_1A3_1DNA", "CoM_01T_1A3_1DNA", "HT29_01Z_1A3_1DNA",
]

P3_SAMPLES = P1_SAMPLES.copy()  # Same 16 as P1

P4_SAMPLES = [
    "CoB_02W_1A3_1RNA", "CoB_08R_3A2_TNA-mRT-EM", "CoB_08X_3A2_TNA-RT-EM",
    "CoB_09D_3A2_RNA-mRT-EM", "CoB_09J_3A2_RNA-RT-EM",
    "CoM_02V_1A3_1RNA", "CoM_08S_3A2_TNA-mRT-EM", "CoM_08Y_3A2_TNA-RT-EM",
    "CoM_09E_3A2_RNA-mRT-EM", "CoM_09K_3A2_RNA-RT-EM",
    "HT29_02T_1A3_1RNA",
    "HT29_21S_3A2_bsTNA-HP", "HT29_21T_3A2_bsTNA-HP",
    "HT29_21U_3A2_bsRNA-HP", "HT29_21V_3A2_bsRNA-HP",
]

P5_SAMPLES = P4_SAMPLES.copy()  # Same 15 as P4

# Key output files per pipeline (file suffix patterns)
# Format: (column_name, file_suffix)
P1_FILES = [
    ("vcf_gz", ".bcftools.vcf.gz"),
    ("vcf_gz_tbi", ".bcftools.vcf.gz.tbi"),
    ("revelio_bam", ".revelio.sorted.bam"),
    ("revelio_bam_bai", ".revelio.sorted.bam.bai"),
    ("sorted_bam", ".sorted.bam"),
    ("markdup_metrics", ".markdup_metrics.txt"),
]

P2_FILES = [
    ("methylation_vcf_gz", ".methylation.vcf.gz"),
    ("methylation_vcf_gz_tbi", ".methylation.vcf.gz.tbi"),
    ("markdup_bam", ".markdup.bam"),
    ("markdup_bam_bai", ".markdup.bam.bai"),
    ("markdup_metrics", ".markdup_metrics.txt"),
]

P3_FILES = [
    ("cna_seg", ".cna.seg"),
    ("seg", ".seg"),
    ("seg_txt", ".seg.txt"),
    ("correctedDepth_txt", ".correctedDepth.txt"),
    ("params_txt", ".params.txt"),
    ("RData", ".RData"),
    ("wig", ".wig"),
]

P4_FILES = [
    ("featureCounts_tsv", ".featureCounts.tsv"),
    ("featureCounts_summary", ".featureCounts.tsv.summary"),
    ("aligned_bam", ".Aligned.sortedByCoord.out.bam"),
]

P5_FILES = [
    ("vcf_gz", ".bcftools.vcf.gz"),
    ("vcf_gz_tbi", ".bcftools.vcf.gz.tbi"),
    ("revelio_bam", ".revelio.bam"),
    ("revelio_bam_bai", ".revelio.bam.bai"),
]

# The "primary" file that determines completion status
P1_PRIMARY = ".bcftools.vcf.gz"
P2_PRIMARY = ".methylation.vcf.gz"
P3_PRIMARY = ".cna.seg"
P4_PRIMARY = ".featureCounts.tsv"
P5_PRIMARY = ".bcftools.vcf.gz"

# ── Helper functions ───────────────────────────────────────────────────────────

def s3_ls_recursive(s3_path):
    """List all files under an S3 path recursively. Returns set of full S3 keys."""
    cmd = ["aws", "s3", "ls", "--recursive", s3_path, "--region", REGION]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return set()
        files = set()
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            # Format: "2026-02-10 07:23:31  472084592 Laboratory/SINGLE_V_..."
            parts = line.split(None, 3)
            if len(parts) >= 4:
                files.add("s3://motleybio/" + parts[3])
        return files
    except subprocess.TimeoutExpired:
        print(f"  WARNING: Timeout listing {s3_path}", file=sys.stderr)
        return set()


def s3_ls_top_level(s3_path):
    """List top-level entries (files and prefixes) under an S3 path. Returns set of names."""
    cmd = ["aws", "s3", "ls", s3_path, "--region", REGION]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            return set()
        entries = set()
        for line in result.stdout.strip().split("\n"):
            line = line.strip()
            if not line:
                continue
            if line.startswith("PRE "):
                # Directory prefix
                entries.add(line.replace("PRE ", "").rstrip("/"))
            else:
                parts = line.split(None, 3)
                if len(parts) >= 4:
                    entries.add(parts[3])
        return entries
    except subprocess.TimeoutExpired:
        return set()


def load_sample_metadata(manifest_path):
    """
    Load sample metadata from sample_manifest.csv.
    Returns dict: base_sample_id -> {cell_line, library_type, assay_category, conversion_method, analyte}

    Handles multi-lane samples (e.g. CoB_08L_3A2_DNA-EM_L001) by stripping _L00X suffix
    to create a unique base sample ID.
    """
    metadata = {}
    with open(manifest_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row["sample_id"]
            # Strip lane suffix for multi-lane samples
            base_id = sample_id
            for suffix in ["_L001", "_L002", "_L003", "_L004"]:
                if base_id.endswith(suffix):
                    base_id = base_id[: -len(suffix)]
                    break

            # Only store first occurrence (lanes have same metadata)
            if base_id not in metadata:
                metadata[base_id] = {
                    "cell_line": row["cell_line"],
                    "library_type": row["library_type"],
                    "assay_category": row["assay_category"],
                    "conversion_method": row["conversion_method"],
                    "analyte": row["analyte"],
                }
    return metadata


def determine_status(existing_files, sample, pipeline_dir, primary_suffix):
    """
    Determine sample status:
    - 'complete' if primary output file exists
    - 'running' if sample directory exists but primary file is missing
    - 'pending' if sample directory does not exist at all
    """
    primary_path = f"{S3_BASE}/{pipeline_dir}/{sample}/{sample}{primary_suffix}"
    if primary_path in existing_files:
        return "complete"

    # Check if ANY file exists for this sample (means job started)
    sample_prefix = f"{S3_BASE}/{pipeline_dir}/{sample}/"
    for f in existing_files:
        if f.startswith(sample_prefix):
            return "running"

    return "pending"


def generate_pipeline_manifest(pipeline_name, pipeline_dir, samples, file_specs,
                                primary_suffix, metadata, output_dir):
    """Generate a manifest CSV for one pipeline."""
    print(f"\n{'='*70}")
    print(f"Processing {pipeline_name} ({pipeline_dir})")
    print(f"  Samples: {len(samples)}")
    print(f"  Listing S3 files...")

    # List all existing directories at the pipeline level first
    existing_dirs = s3_ls_top_level(f"{S3_BASE}/{pipeline_dir}/")
    print(f"  Found {len(existing_dirs)} sample directories in S3")

    # For each sample that has a directory, list files recursively
    all_existing_files = set()
    samples_with_dirs = []
    for sample in samples:
        if sample in existing_dirs:
            samples_with_dirs.append(sample)

    print(f"  Listing files for {len(samples_with_dirs)} matching sample directories...")
    for sample in samples_with_dirs:
        files = s3_ls_recursive(f"{S3_BASE}/{pipeline_dir}/{sample}/")
        all_existing_files.update(files)
        sys.stdout.write(f"    {sample}: {len(files)} files\n")
        sys.stdout.flush()
    print(f"  Found {len(all_existing_files)} total files")

    # Build CSV rows
    meta_columns = ["sample", "cell_line", "library_type", "assay_category",
                     "conversion_method", "analyte"]
    file_columns = [spec[0] for spec in file_specs]
    header = meta_columns + file_columns + ["status"]

    rows = []
    status_counts = defaultdict(int)

    for sample in samples:
        meta = metadata.get(sample, {})
        row = {
            "sample": sample,
            "cell_line": meta.get("cell_line", ""),
            "library_type": meta.get("library_type", ""),
            "assay_category": meta.get("assay_category", ""),
            "conversion_method": meta.get("conversion_method", ""),
            "analyte": meta.get("analyte", ""),
        }

        # Fill in file paths (expected path, empty if not found)
        for col_name, suffix in file_specs:
            expected_path = f"{S3_BASE}/{pipeline_dir}/{sample}/{sample}{suffix}"
            if expected_path in all_existing_files:
                row[col_name] = expected_path
            else:
                row[col_name] = ""

        # Determine status
        status = determine_status(all_existing_files, sample, pipeline_dir, primary_suffix)
        row["status"] = status
        status_counts[status] += 1
        rows.append(row)

    # Write CSV
    output_path = os.path.join(output_dir, f"{pipeline_dir}_manifest.csv")
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)

    print(f"  Wrote: {output_path}")
    print(f"  Status summary:")
    for status in ["complete", "running", "pending"]:
        count = status_counts.get(status, 0)
        print(f"    {status}: {count}/{len(samples)}")

    return output_path, rows, status_counts


def print_final_summary(results):
    """Print a consolidated summary table."""
    print(f"\n{'='*70}")
    print("FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"{'Pipeline':<25} {'Total':>6} {'Complete':>9} {'Running':>8} {'Pending':>8}")
    print(f"{'-'*25} {'-'*6} {'-'*9} {'-'*8} {'-'*8}")
    for name, (path, rows, status_counts) in results.items():
        total = len(rows)
        complete = status_counts.get("complete", 0)
        running = status_counts.get("running", 0)
        pending = status_counts.get("pending", 0)
        print(f"{name:<25} {total:>6} {complete:>9} {running:>8} {pending:>8}")
    print(f"\nOutput files:")
    for name, (path, rows, status_counts) in results.items():
        print(f"  {path}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading sample metadata from manifest...")
    metadata = load_sample_metadata(MANIFEST_PATH)
    print(f"  Loaded metadata for {len(metadata)} unique samples")

    # Pipeline definitions
    pipelines = [
        ("P1 DNA SNP",    "p1_dna_snp",    P1_SAMPLES, P1_FILES, P1_PRIMARY),
        ("P2 DNA Meth",   "p2_dna_meth",   P2_SAMPLES, P2_FILES, P2_PRIMARY),
        ("P3 CNV",        "p3_cnv",        P3_SAMPLES, P3_FILES, P3_PRIMARY),
        ("P4 RNA Counts", "p4_rna_counts", P4_SAMPLES, P4_FILES, P4_PRIMARY),
        ("P5 RNA SNP",    "p5_rna_snp",    P5_SAMPLES, P5_FILES, P5_PRIMARY),
    ]

    results = {}
    for name, pipeline_dir, samples, file_specs, primary_suffix in pipelines:
        path, rows, status_counts = generate_pipeline_manifest(
            name, pipeline_dir, samples, file_specs, primary_suffix, metadata, OUTPUT_DIR
        )
        results[name] = (path, rows, status_counts)

    print_final_summary(results)

    # Print detailed per-sample status for incomplete pipelines
    print(f"\n{'='*70}")
    print("INCOMPLETE SAMPLES DETAIL")
    print(f"{'='*70}")
    any_incomplete = False
    for name, (path, rows, status_counts) in results.items():
        incomplete = [r for r in rows if r["status"] != "complete"]
        if incomplete:
            any_incomplete = True
            print(f"\n{name}:")
            for r in incomplete:
                print(f"  {r['sample']}: {r['status']}")
    if not any_incomplete:
        print("  All samples complete across all pipelines!")


if __name__ == "__main__":
    main()
