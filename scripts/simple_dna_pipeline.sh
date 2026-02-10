#!/bin/bash
# Simple DNA Pipeline - P1 (SNP), P2 (Methylation), P3 (CNV)
# Shared alignment with Biscuit + MarkDuplicates, then three outputs
# No Nextflow, just straightforward bash
#
# Key changes:
# - BCFtools replaces LoFreq (50-100x faster)
# - Conditional Revelio (only for bisulfite/EM treated samples)
# - Downsampling support for SingleAnalyte samples
# - BWA option for non-bisulfite samples

set -e

# Configuration
THREADS=40
WORKDIR=/data/dna_work
RESULTS=/data/results
REF=~/references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa
BISCUIT_INDEX=/data/references/biscuit_index/GRCh38_full_analysis_set_plus_decoy_hla.fa
S3_OUTPUT="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"
REVELIO_SCRIPT=$(dirname "$0")/../pipelines/modules/local/revelio/bin/revelio.py
REVELIO_PARALLEL=$(dirname "$0")/revelio_parallel.sh
DOWNSAMPLE_SCRIPT=$(dirname "$0")/downsample_bam.sh

# Downsampling targets
DOWNSAMPLE_DNA=140000000  # 140M reads for SingleAnalyte DNA to match mTNA

# Container images
BISCUIT_IMG="community.wave.seqera.io/library/biscuit_samtools:84373c8a97fa63b8"
# PICARD_IMG removed - using samtools markdup instead
CNVNATOR_IMG="quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11"

# Create directories
mkdir -p "$WORKDIR" "$RESULTS/p1_dna_snp" "$RESULTS/p2_dna_meth" "$RESULTS/p3_cnv"

echo "=========================================="
echo "DNA Pipeline (P1 + P2 + P3) - BCFtools"
echo "=========================================="
echo "Threads: $THREADS"
echo ""

# Check biscuit index exists
check_biscuit_index() {
    if [ ! -f "$BISCUIT_INDEX" ]; then
        echo "Biscuit index not found. Downloading from S3..."
        mkdir -p /data/references/biscuit_index
        aws s3 cp --quiet s3://motleybio/Resources/biscuit_reference_genome/ /data/references/biscuit_index/ --recursive
    fi
}

# All samples get Revelio for pipeline consistency
# Even non-bisulfite samples (WGS) - ensures identical processing so differences reflect biology
needs_revelio() {
    # Always return true - apply Revelio to ALL samples for comparability
    return 0
}

# Determine if sample needs downsampling
needs_downsampling() {
    local SAMPLE=$1
    case "$SAMPLE" in
        *"1A3"*)
            # SingleAnalyte samples (1A3 in name) need downsampling
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

# Determine if sample should use BWA instead of Biscuit
# BWA is better for WGS (non-bisulfite) samples
use_bwa() {
    local SAMPLE=$1
    case "$SAMPLE" in
        *"WGS"*)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

# Process a single DNA sample through the full pipeline
process_dna_sample() {
    local SAMPLE=$1
    local FQ1=$2
    local FQ2=$3  # Empty for single-end
    local SINGLE_END=$4

    echo ""
    echo "=========================================="
    echo "Processing: $SAMPLE"
    echo "=========================================="

    local SAMPLE_DIR="$WORKDIR/$SAMPLE"
    mkdir -p "$SAMPLE_DIR"

    local ALIGNED_BAM="$SAMPLE_DIR/${SAMPLE}.aligned.bam"
    local SORTED_BAM="$SAMPLE_DIR/${SAMPLE}.sorted.bam"
    local MARKDUP_BAM="$SAMPLE_DIR/${SAMPLE}.markdup.bam"
    local CALMD_BAM="$SAMPLE_DIR/${SAMPLE}.calmd.bam"

    # Step 0: Downsample FASTQs if needed (for SingleAnalyte samples)
    local WORKING_FQ1="$FQ1"
    local WORKING_FQ2="$FQ2"
    if needs_downsampling "$SAMPLE"; then
        echo "[0/4] SingleAnalyte sample - downsampling will be done post-alignment"
        # Note: For DNA, we downsample post-alignment to avoid seqtk dependency
        # and to get more accurate read counts
    fi

    # Step 1: Alignment (Biscuit for bisulfite, BWA for WGS)
    if [ ! -f "$SORTED_BAM" ]; then
        if use_bwa "$SAMPLE"; then
            echo "[1/4] Running BWA-MEM2 alignment (WGS sample)..."
            if [ "$SINGLE_END" = "true" ]; then
                bwa-mem2 mem -t $THREADS "$REF" "$WORKING_FQ1" \
                    | samtools sort -@ 8 -o "$SORTED_BAM" -
            else
                bwa-mem2 mem -t $THREADS "$REF" "$WORKING_FQ1" "$WORKING_FQ2" \
                    | samtools sort -@ 8 -o "$SORTED_BAM" -
            fi
        else
            echo "[1/4] Running Biscuit alignment (bisulfite/EM sample)..."
            if [ "$SINGLE_END" = "true" ]; then
                docker run --rm \
                    -v "$WORKDIR:$WORKDIR" \
                    -v "$(dirname $BISCUIT_INDEX):$(dirname $BISCUIT_INDEX)" \
                    -v "$(dirname $FQ1):$(dirname $FQ1)" \
                    $BISCUIT_IMG \
                    bash -c "biscuit align -@ $THREADS $BISCUIT_INDEX $WORKING_FQ1 | samtools sort -@ 8 -o $SORTED_BAM -"
            else
                docker run --rm \
                    -v "$WORKDIR:$WORKDIR" \
                    -v "$(dirname $BISCUIT_INDEX):$(dirname $BISCUIT_INDEX)" \
                    -v "$(dirname $FQ1):$(dirname $FQ1)" \
                    $BISCUIT_IMG \
                    bash -c "biscuit align -@ $THREADS $BISCUIT_INDEX $WORKING_FQ1 $WORKING_FQ2 | samtools sort -@ 8 -o $SORTED_BAM -"
            fi
        fi
        samtools index -@ $THREADS "$SORTED_BAM"
    else
        echo "[1/4] Aligned BAM exists, skipping..."
    fi

    # Step 1b: Downsample aligned BAM if needed (SingleAnalyte)
    local WORKING_BAM="$SORTED_BAM"
    if needs_downsampling "$SAMPLE"; then
        local DOWNSAMPLED_BAM="$SAMPLE_DIR/${SAMPLE}.downsampled.bam"
        if [ ! -f "$DOWNSAMPLED_BAM" ]; then
            echo "[1b/4] Downsampling to ${DOWNSAMPLE_DNA} reads..."
            "$DOWNSAMPLE_SCRIPT" -i "$SORTED_BAM" -o "$DOWNSAMPLED_BAM" -t "$DOWNSAMPLE_DNA" -@ "$THREADS"
        else
            echo "[1b/4] Downsampled BAM exists, skipping..."
        fi
        WORKING_BAM="$DOWNSAMPLED_BAM"
    fi

    # Step 2: Mark duplicates with samtools markdup
    if [ ! -f "$MARKDUP_BAM" ]; then
        echo "[2/4] Running samtools markdup..."
        samtools collate -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$WORKING_BAM"
        samtools fixmate -@ $THREADS -m "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
        samtools sort -@ $THREADS -o "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam"
        samtools markdup -@ $THREADS -s "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam" "$MARKDUP_BAM" 2> "$SAMPLE_DIR/${SAMPLE}.markdup_metrics.txt"
        samtools index -@ $THREADS "$MARKDUP_BAM"
        # Clean up intermediates
        rm -f "$SAMPLE_DIR/${SAMPLE}.collate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.bam" "$SAMPLE_DIR/${SAMPLE}.fixmate.sorted.bam"
    else
        echo "[2/4] MarkDup BAM exists, skipping..."
    fi

    # Step 3: Add MD tags with samtools calmd
    if [ ! -f "$CALMD_BAM" ]; then
        echo "[3/4] Running samtools calmd..."
        samtools calmd -b -@ $THREADS "$MARKDUP_BAM" "$REF" > "$CALMD_BAM" 2> "$SAMPLE_DIR/calmd.log"
        samtools index -@ $THREADS "$CALMD_BAM"
    else
        echo "[3/4] calmd BAM exists, skipping..."
    fi

    # Step 4a: P1 - Revelio (conditional) + BCFtools for SNP calling
    run_p1_snp "$SAMPLE" "$CALMD_BAM"

    # Step 4b: P2 - Biscuit pileup for methylation (only for bisulfite samples)
    if ! use_bwa "$SAMPLE"; then
        run_p2_meth "$SAMPLE" "$MARKDUP_BAM"
    else
        echo "  [P2] Skipping methylation for WGS sample"
    fi

    # Step 4c: P3 - CNVpytor for copy number
    run_p3_cnv "$SAMPLE" "$MARKDUP_BAM"

    echo "Done: $SAMPLE"
}

run_p1_snp() {
    local SAMPLE=$1
    local CALMD_BAM=$2
    local OUT_DIR="$RESULTS/p1_dna_snp/$SAMPLE"
    mkdir -p "$OUT_DIR"

    local REVELIO_BAM="$WORKDIR/$SAMPLE/${SAMPLE}.revelio.bam"
    local REVELIO_SORTED="$WORKDIR/$SAMPLE/${SAMPLE}.revelio.sorted.bam"
    local VCF="$OUT_DIR/${SAMPLE}.bcftools.vcf"

    if [ -f "$VCF" ]; then
        echo "  [P1] VCF exists, skipping..."
        return
    fi

    local FINAL_BAM="$REVELIO_SORTED"

    # Run Revelio on ALL samples for pipeline consistency (PARALLEL)
    if [ ! -f "$REVELIO_SORTED" ]; then
        echo "  [P1] Running Parallel Revelio (12 chunks, applied to ALL samples)..."
        "$REVELIO_PARALLEL" -i "$CALMD_BAM" -o "$REVELIO_SORTED" -n 12 -t 3
    else
        echo "  [P1] Revelio BAM exists, skipping..."
    fi

    # Run BCFtools for variant calling
    # DNA SNP uses higher quality thresholds: -q 20 -Q 20
    echo "  [P1] Running BCFtools mpileup + call..."
    bcftools mpileup \
        --threads $THREADS \
        -Ou \
        -q 20 \
        -Q 20 \
        -f "$REF" \
        "$FINAL_BAM" \
    | bcftools call \
        --threads $THREADS \
        -mv \
        -Ov \
        -o "$VCF"

    # Compress and index VCF
    bgzip -c "$VCF" > "${VCF}.gz"
    tabix -p vcf "${VCF}.gz"

    echo "  [P1] Done: $VCF"
}

run_p2_meth() {
    local SAMPLE=$1
    local MARKDUP_BAM=$2
    local OUT_DIR="$RESULTS/p2_dna_meth/$SAMPLE"
    mkdir -p "$OUT_DIR"

    local METH_VCF="$OUT_DIR/${SAMPLE}.methylation.vcf"
    local METH_BED="$OUT_DIR/${SAMPLE}.methylation.bed"

    if [ -f "$METH_VCF" ]; then
        echo "  [P2] Methylation VCF exists, skipping..."
        return
    fi

    echo "  [P2] Running Biscuit pileup..."
    docker run --rm \
        -v "$WORKDIR:$WORKDIR" \
        -v "$RESULTS:$RESULTS" \
        -v "$(dirname $BISCUIT_INDEX):$(dirname $BISCUIT_INDEX)" \
        $BISCUIT_IMG \
        bash -c "biscuit pileup -@ $THREADS $BISCUIT_INDEX $MARKDUP_BAM -o $METH_VCF"

    # Also create BED format for easier analysis
    docker run --rm \
        -v "$RESULTS:$RESULTS" \
        -v "$(dirname $BISCUIT_INDEX):$(dirname $BISCUIT_INDEX)" \
        $BISCUIT_IMG \
        bash -c "biscuit vcf2bed -t cg $METH_VCF > $METH_BED"

    echo "  [P2] Done: $METH_VCF"
}

run_p3_cnv() {
    local SAMPLE=$1
    local MARKDUP_BAM=$2
    local OUT_DIR="$RESULTS/p3_cnv/$SAMPLE"
    mkdir -p "$OUT_DIR"

    local CNV_ROOT="$OUT_DIR/${SAMPLE}.cnvnator.root"
    local CNV_CALLS="$OUT_DIR/${SAMPLE}.cnv_calls.txt"

    if [ -f "$CNV_CALLS" ]; then
        echo "  [P3] CNV calls exist, skipping..."
        return
    fi

    echo "  [P3] Running CNVnator..."
    # CNVnator requires a specific workflow
    docker run --rm \
        -v "$WORKDIR:$WORKDIR" \
        -v "$RESULTS:$RESULTS" \
        -v "$(dirname $REF):$(dirname $REF)" \
        $CNVNATOR_IMG \
        bash -c "
            cd $OUT_DIR
            cnvnator -root $CNV_ROOT -tree $MARKDUP_BAM
            cnvnator -root $CNV_ROOT -his 1000 -fasta $REF
            cnvnator -root $CNV_ROOT -stat 1000
            cnvnator -root $CNV_ROOT -partition 1000
            cnvnator -root $CNV_ROOT -call 1000 > $CNV_CALLS
        "

    echo "  [P3] Done: $CNV_CALLS"
}

# Download sample if needed and return local path
download_sample() {
    local SAMPLE=$1
    local S3_PATH=$2
    local LOCAL_PATH="/data/fastq/$(basename $S3_PATH)"

    if [ ! -f "$LOCAL_PATH" ]; then
        echo "Downloading $S3_PATH..."
        mkdir -p /data/fastq
        aws s3 cp --quiet "$S3_PATH" "$LOCAL_PATH"
    fi
    echo "$LOCAL_PATH"
}

# TrinitySeq DNA Samples (10 total)
# CoB: DNA-EM, TNA-mRT-EM, TNA-RT-EM (3)
# CoM: DNA-EM, TNA-mRT-EM, TNA-RT-EM (3)
# HT29: bsTNA-HP x2, bsDNA-HP x2 (4)
TRINITYSEQ_DNA_SAMPLES=(
    "CoB_08L_3A2_DNA-EM"
    "CoB_08R_3A2_TNA-mRT-EM"
    "CoB_08X_3A2_TNA-RT-EM"
    "CoM_08M_3A2_DNA-EM"
    "CoM_08S_3A2_TNA-mRT-EM"
    "CoM_08Y_3A2_TNA-RT-EM"
    "HT29_21S_3A2_bsTNA-HP"
    "HT29_21T_3A2_bsTNA-HP"
    "HT29_21W_3A2_bsDNA-HP"
    "HT29_21X_3A2_bsDNA-HP"
)

# SingleAnalyte DNA Samples (6 total) - will be downsampled
SINGLEANALYTE_DNA_SAMPLES=(
    "CoB_02Y_1A3_1WGS"
    "CoB_02X_1A3_1WGEM"
    "CoM_02U_1A3_1WGS"
    "CoM_02Z_1A3_1WGEM"
    "HT29_02R_1A3_1WGS"
    "HT29_02S_1A3_1WGEM"
)

# Main execution
main() {
    check_biscuit_index

    echo ""
    echo "Processing DNA samples..."
    echo ""
    echo "Available sample arrays:"
    echo "  TRINITYSEQ_DNA_SAMPLES: ${#TRINITYSEQ_DNA_SAMPLES[@]} samples"
    echo "  SINGLEANALYTE_DNA_SAMPLES: ${#SINGLEANALYTE_DNA_SAMPLES[@]} samples"
    echo ""
    echo "To process a sample, call:"
    echo "  process_dna_sample \"SAMPLE_NAME\" \"/path/to/R1.fq\" \"/path/to/R2.fq\" \"false\""
    echo ""
    echo "For single-end (HairyTNA):"
    echo "  process_dna_sample \"HT29_21S_3A2_bsTNA-HP\" \"/data/fastq/HT29.fq\" \"\" \"true\""
    echo ""

    # Example: Process one sample
    # process_dna_sample "SAMPLE_NAME" "/path/to/R1.fq" "/path/to/R2.fq" "false"

    # For HairyTNA (single-end):
    # process_dna_sample "HT29_21S_3A2_bsTNA-HP" "/data/fastq/HT29.fq" "" "true"
}

# Allow sourcing this script for function access, or running directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
