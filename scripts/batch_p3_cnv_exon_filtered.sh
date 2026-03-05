#!/bin/bash
# P3 CNV pipeline with exonic region filtering for AWS Batch
#
# Same as batch_p3_cnv.sh but filters out reads in exonic regions
# BEFORE readCounter, so IchorCNA bins only count non-exonic reads.
# This removes potential RNA bleedthrough artifacts in TrinitySeq samples.
#
# Usage: batch_p3_cnv_exon_filtered.sh <SAMPLE> <BAM_PATH> [BAM_TYPE]
#   SAMPLE:    sample name
#   BAM_PATH:  S3 path to BAM file
#   BAM_TYPE:  "markdup" (default) or "sorted" (will run markdup first)
set -euo pipefail

SAMPLE="${1:?Usage: $0 SAMPLE BAM_PATH [BAM_TYPE]}"
BAM_PATH="${2:?Missing BAM_PATH}"
BAM_TYPE="${3:-markdup}"

THREADS=$(nproc)
WORKDIR=/scratch
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
S3_REFS="s3://motleybio/Resources"
REF="$WORKDIR/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Non-exonic BED for read filtering
S3_NON_EXONIC_BED="$S3_OUTPUT/references/non_exonic_regions.hg38.bed"
NON_EXONIC_BED="$WORKDIR/refs/non_exonic_regions.hg38.bed"

# IchorCNA settings — 1Mb bins
BIN_SIZE=1000000

# Locate IchorCNA package paths
ICHORCNA_PKG=$(Rscript --vanilla -e "cat(system.file(package='ichorCNA'))")
EXTDATA="$ICHORCNA_PKG/extdata"

echo "=========================================="
echo "AWS Batch P3 CNV Pipeline (Exon-Filtered)"
echo "=========================================="
echo "Sample: $SAMPLE"
echo "BAM path: $BAM_PATH"
echo "BAM type: $BAM_TYPE"
echo "Bin size: $((BIN_SIZE / 1000))kb"
echo "Threads: $THREADS"
echo "IchorCNA pkg: $ICHORCNA_PKG"
echo "Started: $(date)"
echo ""

# ---- Step 0: Clean stale state from previous jobs ----
rm -rf "$WORKDIR/bams" "$WORKDIR/results"

# ---- Step 1: Download references ----
echo "[1/5] Downloading references..."
mkdir -p "$WORKDIR/refs" "$WORKDIR/bams" "$WORKDIR/results"

if [ ! -f "$REF" ]; then
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" "$REF"
    aws s3 cp --quiet "$S3_REFS/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" "${REF}.fai"
fi

if [ ! -f "$NON_EXONIC_BED" ]; then
    aws s3 cp --quiet "$S3_NON_EXONIC_BED" "$NON_EXONIC_BED"
fi
echo "References ready."

# ---- Step 2: Get analysis-ready BAM ----
ANALYSIS_BAM="$WORKDIR/bams/${SAMPLE}.analysis.bam"

echo "[2/5] Downloading BAM..."
DOWNLOADED_BAM="$WORKDIR/bams/${SAMPLE}.downloaded.bam"
aws s3 cp --quiet "$BAM_PATH" "$DOWNLOADED_BAM"

# Download or create BAI
BAI_PATH="${BAM_PATH}.bai"
if aws s3 ls "$BAI_PATH" > /dev/null 2>&1; then
    aws s3 cp --quiet "$BAI_PATH" "${DOWNLOADED_BAM}.bai"
else
    echo "  No BAI found, indexing..."
    samtools index -@ $THREADS "$DOWNLOADED_BAM"
fi

if [ "$BAM_TYPE" = "sorted" ]; then
    echo "  Running markdup on sorted BAM..."
    samtools collate -@ $THREADS -o "$WORKDIR/bams/collate.bam" "$DOWNLOADED_BAM"
    samtools fixmate -@ $THREADS -m "$WORKDIR/bams/collate.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/collate.bam"
    samtools sort -@ $THREADS -o "$WORKDIR/bams/fixmate.sorted.bam" "$WORKDIR/bams/fixmate.bam"
    rm -f "$WORKDIR/bams/fixmate.bam"
    samtools markdup -@ $THREADS -s "$WORKDIR/bams/fixmate.sorted.bam" "$ANALYSIS_BAM" 2> "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt"
    rm -f "$WORKDIR/bams/fixmate.sorted.bam"
    samtools index -@ $THREADS "$ANALYSIS_BAM"
    rm -f "$DOWNLOADED_BAM" "${DOWNLOADED_BAM}.bai"
    echo "  Markdup done."
else
    mv "$DOWNLOADED_BAM" "$ANALYSIS_BAM"
    mv "${DOWNLOADED_BAM}.bai" "${ANALYSIS_BAM}.bai" 2>/dev/null || samtools index -@ $THREADS "$ANALYSIS_BAM"
    echo "  Using markdup BAM directly."
fi

echo "Analysis BAM ready: $(ls -lh $ANALYSIS_BAM)"

# ---- Step 3: Filter out exonic reads ----
FILTERED_BAM="$WORKDIR/bams/${SAMPLE}.non_exonic.bam"

echo "[3/5] Filtering out exonic reads..."
TOTAL_READS=$(samtools view -c "$ANALYSIS_BAM")
echo "  Total reads before filtering: $TOTAL_READS"

samtools view -b -L "$NON_EXONIC_BED" -@ $THREADS "$ANALYSIS_BAM" > "$FILTERED_BAM"
samtools index -@ $THREADS "$FILTERED_BAM"

FILTERED_READS=$(samtools view -c "$FILTERED_BAM")
REMOVED_READS=$((TOTAL_READS - FILTERED_READS))
if [ "$TOTAL_READS" -gt 0 ]; then
    PCT_REMOVED=$(awk "BEGIN { printf \"%.1f\", $REMOVED_READS * 100 / $TOTAL_READS }")
else
    PCT_REMOVED="0.0"
fi
echo "  Reads after filtering: $FILTERED_READS (removed $REMOVED_READS = ${PCT_REMOVED}% exonic)"

# Free unfiltered BAM
rm -f "$ANALYSIS_BAM" "${ANALYSIS_BAM}.bai"

# ---- Step 4: IchorCNA ----
OUTDIR="$WORKDIR/results/ichorcna"
mkdir -p "$OUTDIR"

echo "[4/5] Running IchorCNA pipeline on exon-filtered BAM..."

# Step 4a: readCounter
echo "  Running readCounter (${BIN_SIZE}bp bins)..."
CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
readCounter \
    --window "$BIN_SIZE" \
    --quality 20 \
    --chromosome "$CHRS" \
    "$FILTERED_BAM" > "$OUTDIR/${SAMPLE}.wig"

WIG_LINES=$(wc -l < "$OUTDIR/${SAMPLE}.wig")
echo "  Read counts done: $WIG_LINES lines in WIG"

# Free filtered BAM
rm -f "$FILTERED_BAM" "${FILTERED_BAM}.bai"

# Step 4b: Verify IchorCNA reference files
GC_WIG="$EXTDATA/gc_hg38_1000kb.wig"
MAP_WIG="$EXTDATA/map_hg38_1000kb.wig"
CENTROMERE="$EXTDATA/GRCh38.GCA_000001405.2_centromere_acen.txt"
PON="$EXTDATA/HD_ULP_PoN_hg38_1Mb_normAutosomes_median.rds"

echo "  Checking reference files..."
for f in "$GC_WIG" "$MAP_WIG"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Reference file not found: $f"
        ls -la "$EXTDATA/" || echo "  extdata dir not found"
        exit 1
    fi
    echo "    OK: $(basename $f)"
done

CENTROMERE_R="NULL"
if [ -f "$CENTROMERE" ]; then
    CENTROMERE_R="'$CENTROMERE'"
    echo "    OK: $(basename $CENTROMERE)"
fi

PON_R="NULL"
if [ -f "$PON" ]; then
    PON_R="'$PON'"
    echo "    OK: $(basename $PON)"
fi

# Step 4c: Run IchorCNA via R function
echo "  Running IchorCNA HMM..."
Rscript --vanilla -e "
library(ichorCNA)
run_ichorCNA(
    tumor_wig  = '$OUTDIR/${SAMPLE}.wig',
    id         = '$SAMPLE',
    gcWig      = '$GC_WIG',
    mapWig     = '$MAP_WIG',
    centromere = $CENTROMERE_R,
    normal_panel = $PON_R,
    normal     = 'c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)',
    ploidy     = 'c(2, 3)',
    maxCN      = 5,
    includeHOMD      = FALSE,
    estimateNormal   = TRUE,
    estimatePloidy   = TRUE,
    estimateScPrevalence = TRUE,
    scStates   = 'c(1, 3)',
    genomeBuild = 'hg38',
    genomeStyle = 'UCSC',
    chrs       = 'c(1:22, \"X\")',
    chrTrain   = 'c(1:22)',
    chrNormalize = 'c(1:22)',
    cores      = $THREADS,
    outDir     = '$OUTDIR'
)
"

echo "  IchorCNA complete."

# Summarize results
SEGMENTS="$OUTDIR/${SAMPLE}.seg"
PARAMS="$OUTDIR/${SAMPLE}.params.txt"

if [ -f "$PARAMS" ]; then
    echo "  Parameters:"
    head -20 "$PARAMS"
fi

TOTAL_SEGS=0
if [ -f "$SEGMENTS" ]; then
    TOTAL_SEGS=$(tail -n +2 "$SEGMENTS" | wc -l || echo "0")
    echo "  Total segments: $TOTAL_SEGS"
fi

# ---- Step 5: Upload results ----
echo "[5/5] Uploading results to S3..."
aws s3 cp --recursive --quiet "$OUTDIR/" "$S3_OUTPUT/p3_cnv_exon_filtered/${SAMPLE}/"

if [ -f "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt" ]; then
    aws s3 cp --quiet "$WORKDIR/results/${SAMPLE}.markdup_metrics.txt" "$S3_OUTPUT/p3_cnv_exon_filtered/${SAMPLE}/${SAMPLE}.markdup_metrics.txt"
fi

# Write filtering stats for downstream summary
cat > "$WORKDIR/results/${SAMPLE}.filter_stats.txt" <<STATS
sample=$SAMPLE
total_reads=$TOTAL_READS
filtered_reads=$FILTERED_READS
removed_reads=$REMOVED_READS
pct_exonic_removed=$PCT_REMOVED
total_segments=$TOTAL_SEGS
STATS
aws s3 cp --quiet "$WORKDIR/results/${SAMPLE}.filter_stats.txt" "$S3_OUTPUT/p3_cnv_exon_filtered/${SAMPLE}/${SAMPLE}.filter_stats.txt"

echo ""
echo "=========================================="
echo "COMPLETE: $SAMPLE ($TOTAL_SEGS segments)"
echo "  Exonic reads removed: ${PCT_REMOVED}%"
echo "=========================================="
echo "Finished: $(date)"
echo "Results: $S3_OUTPUT/p3_cnv_exon_filtered/${SAMPLE}/"
