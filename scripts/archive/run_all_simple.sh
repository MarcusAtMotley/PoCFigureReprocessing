#!/bin/bash
# Master orchestration script for all pipelines
# Usage:
#   ./run_all_simple.sh setup      # Download references and inputs
#   ./run_all_simple.sh p5         # Run P5 (RNA SNP) only
#   ./run_all_simple.sh dna        # Run DNA pipelines (P1+P2+P3)
#   ./run_all_simple.sh all        # Run everything

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/simple_dna_pipeline.sh"

# DNA Samples configuration
# Format: SAMPLE_NAME|S3_FQ1|S3_FQ2|SINGLE_END
DNA_SAMPLES=(
    # SingleAnalyte WGS (P1 only - no methylation)
    "CoB_02M_1C3_1DNA|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R1.fastq.gz|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_02M_1C3_1DNA_R2.fastq.gz|false"
    "CoM_02K_1C3_1DNA|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R1.fastq.gz|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_02K_1C3_1DNA_R2.fastq.gz|false"
    "HT29_02N_1B3_1DNA|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R1.fastq.gz|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_02N_1B3_1DNA_R2.fastq.gz|false"

    # SingleAnalyte WGEM (P2 methylation)
    "CoB_01W_1A3_1DNA|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R1.fastq.gz|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoB_01W_1A3_1DNA_R2.fastq.gz|false"
    "CoM_01T_1A3_1DNA|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R1.fastq.gz|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/CoM_01T_1A3_1DNA_R2.fastq.gz|false"
    "HT29_01Z_1A3_1DNA|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R1.fastq.gz|s3://motleybio-medgenome/data/P2008501_04082025/FASTQ/HT2_01Z_1A3_1DNA_R2.fastq.gz|false"

    # TrinitySeq DNA-EM (already merged on S3)
    "CoB_08L_3A2_DNA-EM|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoB_08L_3A2_DNA-EM_S1_merged_R1.fastq|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoB_08L_3A2_DNA-EM_S1_merged_R2.fastq|false"
    "CoM_08M_3A2_DNA-EM|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoM_08M_3A2_DNA-EM_S2_merged_R1.fastq|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoM_08M_3A2_DNA-EM_S2_merged_R2.fastq|false"

    # TrinitySeq mTNA (already merged on S3)
    "CoB_08R_3A2_TNA-mRT-EM|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoB_08R_3A2_TNA-mRT-EM_S3_merged_R1.fastq|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoB_08R_3A2_TNA-mRT-EM_S3_merged_R2.fastq|false"
    "CoM_08S_3A2_TNA-mRT-EM|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoM_08S_3A2_TNA-mRT-EM_S4_merged_R1.fastq|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoM_08S_3A2_TNA-mRT-EM_S4_merged_R2.fastq|false"

    # TrinitySeq TNA-RT (already merged on S3)
    "CoB_08X_3A2_TNA-RT-EM|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoB_08X_3A2_TNA-RT-EM_S5_merged_R1.fastq|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoB_08X_3A2_TNA-RT-EM_S5_merged_R2.fastq|false"
    "CoM_08Y_3A2_TNA-RT-EM|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoM_08Y_3A2_TNA-RT-EM_S6_merged_R1.fastq|s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS/merged_fastq/CoM_08Y_3A2_TNA-RT-EM_S6_merged_R2.fastq|false"

    # HairyTNA (single-end)
    "HT29_21S_3A2_bsTNA-HP|s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21S_3A2_bsTNA-HP_S7_unbarcoded.cutadapt.fastq||true"
    "HT29_21T_3A2_bsTNA-HP|s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21T_3A2_bsTNA-HP_S8_unbarcoded.cutadapt.fastq||true"
    "HT29_21W_3A2_bsDNA-HP|s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21W_3A2_bsDNA-HP_S11_unbarcoded.cutadapt.fastq||true"
    "HT29_21X_3A2_bsDNA-HP|s3://motleybio/Laboratory/ASSAY_DEVELOPMENT/run_motley34/seqera_HP_output/rna_deconvolution/cutadapt/HT29_21X_3A2_bsDNA-HP_S12_unbarcoded.cutadapt.fastq||true"
)

setup_references() {
    echo "=========================================="
    echo "Setting up references"
    echo "=========================================="

    # Create /data if it doesn't exist
    sudo mkdir -p /data
    sudo chown ubuntu:ubuntu /data

    mkdir -p /data/references /data/fastq /data/dna_work /data/p5_work /data/results

    # Symlink existing reference
    if [ ! -f /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa ]; then
        if [ -f ~/references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa ]; then
            echo "Symlinking existing reference..."
            ln -sf ~/references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa /data/references/
        else
            echo "Downloading reference FASTA..."
            aws s3 cp --quiet s3://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa /data/references/
        fi
    fi

    # Ensure reference is indexed
    if [ ! -f /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ]; then
        echo "Indexing reference..."
        samtools faidx /data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa
    fi

    # Download biscuit index
    if [ ! -d /data/references/biscuit_index ] || [ -z "$(ls -A /data/references/biscuit_index 2>/dev/null)" ]; then
        echo "Downloading Biscuit index (~10GB)..."
        mkdir -p /data/references/biscuit_index
        aws s3 cp --quiet s3://motleybio/Resources/biscuit_reference_genome/ /data/references/biscuit_index/ --recursive
    else
        echo "Biscuit index already exists"
    fi

    echo "Setup complete!"
    echo ""
    df -h /data
}

download_dna_fastqs() {
    echo "=========================================="
    echo "Downloading DNA FASTQs"
    echo "=========================================="

    mkdir -p /data/fastq

    for entry in "${DNA_SAMPLES[@]}"; do
        IFS='|' read -r sample fq1 fq2 single_end <<< "$entry"
        local_fq1="/data/fastq/$(basename $fq1)"

        if [ ! -f "$local_fq1" ]; then
            echo "Downloading $sample R1..."
            aws s3 cp --quiet "$fq1" "$local_fq1"
        fi

        if [ -n "$fq2" ]; then
            local_fq2="/data/fastq/$(basename $fq2)"
            if [ ! -f "$local_fq2" ]; then
                echo "Downloading $sample R2..."
                aws s3 cp --quiet "$fq2" "$local_fq2"
            fi
        fi
    done

    echo "Download complete!"
    du -sh /data/fastq
}

run_p5() {
    echo "=========================================="
    echo "Running P5 (RNA SNP)"
    echo "=========================================="
    bash "$SCRIPT_DIR/simple_p5_rna_snp.sh"
}

run_dna() {
    echo "=========================================="
    echo "Running DNA Pipelines (P1 + P2 + P3)"
    echo "=========================================="

    # First ensure we have the data
    download_dna_fastqs

    for entry in "${DNA_SAMPLES[@]}"; do
        IFS='|' read -r sample fq1 fq2 single_end <<< "$entry"
        local_fq1="/data/fastq/$(basename $fq1)"
        local_fq2=""
        if [ -n "$fq2" ]; then
            local_fq2="/data/fastq/$(basename $fq2)"
        fi

        process_dna_sample "$sample" "$local_fq1" "$local_fq2" "$single_end"
    done
}

upload_results() {
    echo "=========================================="
    echo "Uploading results to S3"
    echo "=========================================="

    S3_BASE="s3://motleybio/Laboratory/SINGLE_V_TRINITY_COMPARISONS"

    if [ -d /data/results/p1_dna_snp ]; then
        echo "Uploading P1 results..."
        aws s3 sync /data/results/p1_dna_snp/ ${S3_BASE}/p1_dna_snp_results/
    fi

    if [ -d /data/results/p2_dna_meth ]; then
        echo "Uploading P2 results..."
        aws s3 sync /data/results/p2_dna_meth/ ${S3_BASE}/p2_dna_meth_results/
    fi

    if [ -d /data/results/p3_cnv ]; then
        echo "Uploading P3 results..."
        aws s3 sync /data/results/p3_cnv/ ${S3_BASE}/p3_cnv_results/
    fi

    if [ -d /data/results/p5_rna_snp ]; then
        echo "Uploading P5 results..."
        aws s3 sync /data/results/p5_rna_snp/ ${S3_BASE}/p5_rna_snp_results/
    fi

    echo "Upload complete!"
}

show_status() {
    echo "=========================================="
    echo "Current Status"
    echo "=========================================="
    echo ""
    echo "Disk usage:"
    df -h /data 2>/dev/null || df -h /

    echo ""
    echo "Downloaded FASTQs:"
    ls -lh /data/fastq/ 2>/dev/null | head -20 || echo "None yet"

    echo ""
    echo "Results:"
    for dir in p1_dna_snp p2_dna_meth p3_cnv p5_rna_snp; do
        if [ -d "/data/results/$dir" ]; then
            count=$(find /data/results/$dir -name "*.vcf" -o -name "*.txt" 2>/dev/null | wc -l)
            echo "  $dir: $count output files"
        fi
    done
}

usage() {
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  setup     - Download references and create directories"
    echo "  p5        - Run P5 (RNA SNP) pipeline only"
    echo "  dna       - Run DNA pipelines (P1 + P2 + P3)"
    echo "  all       - Run all pipelines"
    echo "  upload    - Upload results to S3"
    echo "  status    - Show current status"
    echo ""
    echo "Recommended order:"
    echo "  1. ./run_all_simple.sh setup"
    echo "  2. ./run_all_simple.sh p5      # Quick win - uses existing P4 BAMs"
    echo "  3. ./run_all_simple.sh dna     # Longer - full alignment + analysis"
    echo "  4. ./run_all_simple.sh upload"
}

case "${1:-}" in
    setup)
        setup_references
        ;;
    p5)
        run_p5
        ;;
    dna)
        run_dna
        ;;
    all)
        setup_references
        run_p5
        run_dna
        ;;
    upload)
        upload_results
        ;;
    status)
        show_status
        ;;
    *)
        usage
        exit 1
        ;;
esac
