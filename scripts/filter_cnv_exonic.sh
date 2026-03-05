#!/bin/bash
# Post-hoc exonic filtering for P3 CNV calls (IchorCNA)
#
# Removes exonic regions from CNV segments to eliminate potential
# RNA bleedthrough artifacts in TrinitySeq samples. Uses GENCODE/Ensembl
# GTF to extract exon coordinates, then bedtools subtract to remove
# exonic base-pairs from each sample's .cna.seg file.
#
# Usage: ./filter_cnv_exonic.sh [--dry-run]
#
# Outputs:
#   - Exon BED:    references/exonic_regions.hg38.bed (merged, sorted)
#   - Per-sample:  <sample>.exon_filtered.cna.seg  (uploaded to S3)
#   - Summary:     p3_cnv_exon_filtered_summary.csv (uploaded to S3)

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="$REPO_DIR/output_manifests/p3_cnv_manifest.csv"
REFS_DIR="$REPO_DIR/references"
WORKDIR="/tmp/cnv_exon_filter"
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/p3_cnv"
S3_GTF="s3://motleybio/Resources/GTF/Homo_sapiens.GRCh38.112.chr_label.gtf"

EXON_BED="$REFS_DIR/exonic_regions.hg38.bed"
DRY_RUN=false

if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "[DRY RUN] Will not upload results to S3"
fi

echo "=========================================="
echo "P3 CNV Exonic Region Filter"
echo "=========================================="
echo "Manifest: $MANIFEST"
echo "Started:  $(date)"
echo ""

mkdir -p "$WORKDIR"

# ---- Step 1: Generate exon BED from GENCODE GTF ----
if [ -f "$EXON_BED" ]; then
    EXON_COUNT=$(wc -l < "$EXON_BED")
    echo "[1/3] Exon BED already exists: $EXON_BED ($EXON_COUNT regions)"
else
    echo "[1/3] Generating exon BED from GENCODE GTF..."
    GTF_LOCAL="$WORKDIR/annotations.gtf"

    if [ ! -f "$GTF_LOCAL" ]; then
        echo "  Downloading GTF from S3..."
        aws s3 cp --quiet "$S3_GTF" "$GTF_LOCAL"
    fi

    # Extract exon features, convert to BED, sort, merge overlapping
    echo "  Extracting exon coordinates..."
    awk '$3 == "exon" { print $1"\t"$4-1"\t"$5 }' "$GTF_LOCAL" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "$EXON_BED"

    EXON_COUNT=$(wc -l < "$EXON_BED")
    echo "  Created $EXON_BED ($EXON_COUNT merged exonic regions)"

    # Clean up GTF (large file)
    rm -f "$GTF_LOCAL"
fi

# Calculate total exonic base-pairs
EXON_BP=$(awk '{ sum += $3 - $2 } END { print sum }' "$EXON_BED")
echo "  Total exonic base-pairs: $(printf "%'d" "$EXON_BP") ($(echo "scale=1; $EXON_BP / 1000000" | bc)Mb)"
echo ""

# ---- Step 2: Filter each sample's CNV calls ----
echo "[2/3] Filtering CNV segments..."
echo ""

# Summary CSV header
SUMMARY_FILE="$WORKDIR/p3_cnv_exon_filtered_summary.csv"
echo "sample,cell_line,assay_category,total_segments_raw,total_bp_raw,total_segments_filtered,total_bp_filtered,bp_removed,pct_bp_removed" > "$SUMMARY_FILE"

# Skip CSV header
SAMPLE_COUNT=0
TOTAL_SAMPLES=$(tail -n +2 "$MANIFEST" | wc -l)

while IFS=',' read -r sample cell_line library_type assay_category conversion_method analyte cna_seg_s3 rest; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo "  [$SAMPLE_COUNT/$TOTAL_SAMPLES] $sample ($assay_category)"

    # Download .cna.seg
    LOCAL_SEG="$WORKDIR/${sample}.cna.seg"
    aws s3 cp --quiet "$cna_seg_s3" "$LOCAL_SEG"

    # Count raw stats (skip header)
    RAW_SEGS=$(tail -n +2 "$LOCAL_SEG" | wc -l)
    RAW_BP=$(tail -n +2 "$LOCAL_SEG" | awk '{ sum += $3 - $2 } END { print sum+0 }')

    # Save header for output
    head -1 "$LOCAL_SEG" > "$WORKDIR/${sample}.header.tmp"

    # Convert seg data to BED-like for bedtools (keep all columns)
    # bedtools subtract needs: chr start end [rest...]
    # .cna.seg already has chr/start/end as first 3 columns
    tail -n +2 "$LOCAL_SEG" > "$WORKDIR/${sample}.body.tmp"

    # bedtools subtract: remove exonic regions from each segment
    # -A flag would remove entire segment if ANY overlap; without it,
    # segments are trimmed/split at exon boundaries
    bedtools subtract \
        -a "$WORKDIR/${sample}.body.tmp" \
        -b "$EXON_BED" \
        > "$WORKDIR/${sample}.filtered.tmp"

    FILT_SEGS=$(wc -l < "$WORKDIR/${sample}.filtered.tmp")
    FILT_BP=$(awk '{ sum += $3 - $2 } END { print sum+0 }' "$WORKDIR/${sample}.filtered.tmp")
    REMOVED_BP=$((RAW_BP - FILT_BP))
    if [ "$RAW_BP" -gt 0 ]; then
        PCT_REMOVED=$(echo "scale=1; $REMOVED_BP * 100 / $RAW_BP" | bc)
    else
        PCT_REMOVED="0.0"
    fi

    # Reassemble with header
    OUTPUT_SEG="$WORKDIR/${sample}.exon_filtered.cna.seg"
    cat "$WORKDIR/${sample}.header.tmp" "$WORKDIR/${sample}.filtered.tmp" > "$OUTPUT_SEG"

    echo "    Raw: $RAW_SEGS segments ($(printf "%'d" "$RAW_BP") bp)"
    echo "    Filtered: $FILT_SEGS segments ($(printf "%'d" "$FILT_BP") bp) — removed ${PCT_REMOVED}% exonic"

    # Upload
    if [ "$DRY_RUN" = false ]; then
        aws s3 cp --quiet "$OUTPUT_SEG" "$S3_OUTPUT/${sample}/${sample}.exon_filtered.cna.seg"
    fi

    # Summary row
    echo "$sample,$cell_line,$assay_category,$RAW_SEGS,$RAW_BP,$FILT_SEGS,$FILT_BP,$REMOVED_BP,$PCT_REMOVED" >> "$SUMMARY_FILE"

    # Clean per-sample temp files
    rm -f "$WORKDIR/${sample}".*.tmp "$LOCAL_SEG"

done < <(tail -n +2 "$MANIFEST")

echo ""

# ---- Step 3: Upload summary ----
echo "[3/3] Summary"
echo ""
column -t -s',' "$SUMMARY_FILE"
echo ""

if [ "$DRY_RUN" = false ]; then
    aws s3 cp --quiet "$SUMMARY_FILE" "$S3_OUTPUT/p3_cnv_exon_filtered_summary.csv"
    # Also keep a local copy
    cp "$SUMMARY_FILE" "$REPO_DIR/output_manifests/p3_cnv_exon_filtered_summary.csv"
    echo "Summary uploaded to: $S3_OUTPUT/p3_cnv_exon_filtered_summary.csv"
fi

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE_COUNT samples filtered"
echo "=========================================="
echo "Exon BED:   $EXON_BED"
echo "Finished:   $(date)"
