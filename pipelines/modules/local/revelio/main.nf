process REVELIO {
    tag "$meta.id"
    label 'process_medium'

    // Wave builds container from conda environment (pysam + samtools)
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.revelio.bam"), emit: bam
    tuple val(meta), path("*.revelio.bam.bai"), emit: bai
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Add MD tags using samtools calmd (required for revelio)
    samtools calmd -b ${bam} ${fasta} > ${prefix}.calmd.bam 2>/dev/null

    # Run revelio to mask bisulfite conversions
    revelio.py ${prefix}.calmd.bam ${prefix}.revelio.bam

    # Index the output BAM
    samtools index ${prefix}.revelio.bam

    # Clean up intermediate files
    rm -f ${prefix}.calmd.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revelio: 1.1.0
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.revelio.bam
    touch ${prefix}.revelio.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revelio: stub
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
