process INTERGENIC_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(cnv_calls)
    path(intergenic_bed)

    output:
    tuple val(meta), path("*_filtered_cnv.csv"), emit: filtered_calls
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Convert CNVpytor TSV output to BED format for intersection
    # CNVpytor output: chr start end size type rd_ratio p_val e_val
    # Skip header if present, convert to BED
    awk 'BEGIN{OFS="\\t"} NR>0 && \$1 ~ /^chr/ {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8}' ${cnv_calls} > cnv_calls.bed

    # Filter CNV calls to intergenic regions (remove gene-overlapping CNVs)
    if [ -s cnv_calls.bed ]; then
        bedtools intersect -a cnv_calls.bed -b ${intergenic_bed} -wa ${args} > filtered.bed

        # Convert back to CSV with header
        echo "chr,start,end,size,type,rd_ratio,p_val,e_val" > ${prefix}_filtered_cnv.csv
        awk 'BEGIN{OFS=","} {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8}' filtered.bed >> ${prefix}_filtered_cnv.csv
    else
        # No CNV calls, create empty output
        echo "chr,start,end,size,type,rd_ratio,p_val,e_val" > ${prefix}_filtered_cnv.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "chr,start,end,size,type,rd_ratio,p_val,e_val" > ${prefix}_filtered_cnv.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
