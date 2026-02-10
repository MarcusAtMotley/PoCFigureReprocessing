#!/bin/bash
# P1 DNA SNP Recovery Script
# After VM crash:
# - Samples 1 & 3 have sorted BAMs, resume from markdup
# - Sample 2 lost its BAM, re-align from FASTQs
# Strategy: Run sample 2 alignment at 16 threads while samples 1&3 do post-alignment at 16 threads each
set -e

SCRIPT_DIR=$(dirname "$0")
source "$SCRIPT_DIR/simple_dna_pipeline.sh"

# Override threads for parallel processing
THREADS=16

echo "=========================================="
echo "P1 DNA SNP - Recovery"
echo "=========================================="
echo "Started: $(date)"
echo ""

# ---- Phase 1: Parallel work ----
# Sample 2: Re-align (longest task, start first)
# Samples 1 & 3: markdup → calmd (can run in parallel)

echo "=== Phase 1: Parallel processing ==="
echo "  Sample 2: Re-aligning (16 threads)"
echo "  Samples 1 & 3: markdup + calmd (16 threads each)"
echo ""

# Sample 2: Re-alignment
(
    SAMPLE="CoM_08M_3A2_DNA-EM"
    SAMPLE_DIR="$WORKDIR/$SAMPLE"
    SORTED_BAM="$SAMPLE_DIR/${SAMPLE}.sorted.bam"
    FQ1="/data/fastq/CoM_08M_3A2_DNA-EM_merged_R1.fastq"
    FQ2="/data/fastq/CoM_08M_3A2_DNA-EM_merged_R2.fastq"

    mkdir -p "$SAMPLE_DIR"
    echo "[$SAMPLE] Starting Biscuit alignment..."
    docker run --rm \
        -v "$WORKDIR:$WORKDIR" \
        -v "$(dirname $BISCUIT_INDEX):$(dirname $BISCUIT_INDEX)" \
        -v "/data/fastq:/data/fastq" \
        $BISCUIT_IMG \
        bash -c "biscuit align -@ $THREADS $BISCUIT_INDEX $FQ1 $FQ2 | samtools sort -@ 8 -o $SORTED_BAM -"
    samtools index -@ $THREADS "$SORTED_BAM"
    echo "[$SAMPLE] Alignment complete!"
) > /data/dna_work/CoM_08M_recovery.log 2>&1 &
ALIGN_PID=$!
echo "Sample 2 alignment started (PID: $ALIGN_PID)"

# Sample 1: markdup → calmd
(
    SAMPLE="CoB_08L_3A2_DNA-EM"
    SAMPLE_DIR="$WORKDIR/$SAMPLE"
    SORTED_BAM="$SAMPLE_DIR/${SAMPLE}.sorted.bam"
    MARKDUP_BAM="$SAMPLE_DIR/${SAMPLE}.markdup.bam"
    CALMD_BAM="$SAMPLE_DIR/${SAMPLE}.calmd.bam"

    if [ ! -f "$MARKDUP_BAM" ]; then
        echo "[$SAMPLE] Running samtools markdup..."
        samtools collate -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SORTED_BAM"
        samtools fixmate -@ $THREADS -m "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
        samtools sort -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
        samtools markdup -@ $THREADS -s "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$SAMPLE_DIR/${SAMPLE}.markdup_metrics.txt"
        samtools index -@ $THREADS "$MARKDUP_BAM"
        rm -f "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam"
        echo "[$SAMPLE] markdup complete"
    else
        echo "[$SAMPLE] markdup BAM exists, skipping"
    fi

    if [ ! -f "$CALMD_BAM" ]; then
        echo "[$SAMPLE] Running calmd..."
        samtools calmd -b -@ $THREADS "$MARKDUP_BAM" ~/references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa > "$CALMD_BAM" 2> "$SAMPLE_DIR/calmd.log"
        samtools index -@ $THREADS "$CALMD_BAM"
        echo "[$SAMPLE] calmd complete"
    else
        echo "[$SAMPLE] calmd BAM exists, skipping"
    fi
) > /data/dna_work/CoB_08L_recovery.log 2>&1 &
S1_PID=$!
echo "Sample 1 markdup+calmd started (PID: $S1_PID)"

# Sample 3: markdup → calmd
(
    SAMPLE="CoB_08R_3A2_TNA-mRT-EM"
    SAMPLE_DIR="$WORKDIR/$SAMPLE"
    SORTED_BAM="$SAMPLE_DIR/${SAMPLE}.sorted.bam"
    MARKDUP_BAM="$SAMPLE_DIR/${SAMPLE}.markdup.bam"
    CALMD_BAM="$SAMPLE_DIR/${SAMPLE}.calmd.bam"

    if [ ! -f "$MARKDUP_BAM" ]; then
        echo "[$SAMPLE] Running samtools markdup..."
        samtools collate -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SORTED_BAM"
        samtools fixmate -@ $THREADS -m "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
        samtools sort -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
        samtools markdup -@ $THREADS -s "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$SAMPLE_DIR/${SAMPLE}.markdup_metrics.txt"
        samtools index -@ $THREADS "$MARKDUP_BAM"
        rm -f "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam"
        echo "[$SAMPLE] markdup complete"
    else
        echo "[$SAMPLE] markdup BAM exists, skipping"
    fi

    if [ ! -f "$CALMD_BAM" ]; then
        echo "[$SAMPLE] Running calmd..."
        samtools calmd -b -@ $THREADS "$MARKDUP_BAM" ~/references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa > "$CALMD_BAM" 2> "$SAMPLE_DIR/calmd.log"
        samtools index -@ $THREADS "$CALMD_BAM"
        echo "[$SAMPLE] calmd complete"
    else
        echo "[$SAMPLE] calmd BAM exists, skipping"
    fi
) > /data/dna_work/CoB_08R_recovery.log 2>&1 &
S3_PID=$!
echo "Sample 3 markdup+calmd started (PID: $S3_PID)"

echo ""
echo "Waiting for samples 1 & 3 markdup+calmd to finish..."
wait $S1_PID && echo "Sample 1 markdup+calmd done" || echo "WARNING: Sample 1 markdup+calmd failed"
wait $S3_PID && echo "Sample 3 markdup+calmd done" || echo "WARNING: Sample 3 markdup+calmd failed"

# ---- Phase 2: Revelio + BCFtools for samples 1 & 3 ----
# Run these while sample 2 is still aligning
echo ""
echo "=== Phase 2: Revelio + BCFtools for samples 1 & 3 ==="

# Sample 1: Revelio + BCFtools (use 20 threads since sample 2 is still aligning)
(
    SAMPLE="CoB_08L_3A2_DNA-EM"
    CALMD_BAM="$WORKDIR/$SAMPLE/${SAMPLE}.calmd.bam"
    echo "[$SAMPLE] Starting Revelio + BCFtools..."
    THREADS=20 run_p1_snp "$SAMPLE" "$CALMD_BAM"
    echo "[$SAMPLE] P1 COMPLETE!"
) > /data/dna_work/CoB_08L_p1.log 2>&1 &
S1_P1_PID=$!
echo "Sample 1 Revelio+BCFtools started"

# Sample 3: Revelio + BCFtools
(
    SAMPLE="CoB_08R_3A2_TNA-mRT-EM"
    CALMD_BAM="$WORKDIR/$SAMPLE/${SAMPLE}.calmd.bam"
    echo "[$SAMPLE] Starting Revelio + BCFtools..."
    THREADS=20 run_p1_snp "$SAMPLE" "$CALMD_BAM"
    echo "[$SAMPLE] P1 COMPLETE!"
) > /data/dna_work/CoB_08R_p1.log 2>&1 &
S3_P1_PID=$!
echo "Sample 3 Revelio+BCFtools started"

echo ""
echo "Waiting for all remaining tasks..."
wait $S1_P1_PID && echo "Sample 1 P1 COMPLETE" || echo "WARNING: Sample 1 P1 failed"
wait $S3_P1_PID && echo "Sample 3 P1 COMPLETE" || echo "WARNING: Sample 3 P1 failed"

# ---- Phase 3: Wait for sample 2 alignment, then process it ----
echo ""
echo "=== Phase 3: Waiting for sample 2 alignment ==="
wait $ALIGN_PID && echo "Sample 2 alignment complete!" || { echo "ERROR: Sample 2 alignment failed!"; exit 1; }

# Now process sample 2 through the full post-alignment pipeline (all 40 threads)
echo ""
echo "=== Phase 4: Processing sample 2 through markdup → calmd → Revelio → BCFtools ==="
THREADS=40
SAMPLE="CoM_08M_3A2_DNA-EM"
SAMPLE_DIR="$WORKDIR/$SAMPLE"
SORTED_BAM="$SAMPLE_DIR/${SAMPLE}.sorted.bam"
MARKDUP_BAM="$SAMPLE_DIR/${SAMPLE}.markdup.bam"
CALMD_BAM="$SAMPLE_DIR/${SAMPLE}.calmd.bam"

echo "[$SAMPLE] Running samtools markdup..."
samtools collate -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SORTED_BAM"
samtools fixmate -@ $THREADS -m "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
samtools sort -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
samtools markdup -@ $THREADS -s "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$SAMPLE_DIR/${SAMPLE}.markdup_metrics.txt"
samtools index -@ $THREADS "$MARKDUP_BAM"
rm -f "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam"
echo "[$SAMPLE] markdup complete"

echo "[$SAMPLE] Running calmd..."
samtools calmd -b -@ $THREADS "$MARKDUP_BAM" ~/references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa > "$CALMD_BAM" 2> "$SAMPLE_DIR/calmd.log"
samtools index -@ $THREADS "$CALMD_BAM"
echo "[$SAMPLE] calmd complete"

echo "[$SAMPLE] Running Revelio + BCFtools..."
run_p1_snp "$SAMPLE" "$CALMD_BAM"
echo "[$SAMPLE] P1 COMPLETE!"

echo ""
echo "=========================================="
echo "P1 Recovery Complete!"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Results:"
ls -la "$RESULTS/p1_dna_snp/"*/
