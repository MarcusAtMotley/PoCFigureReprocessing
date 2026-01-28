process REVELIO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3571571e227f515c2eab-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3571571e227f515c2eab-0' }"

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
    # Clone revelio if not present
    git clone --depth 1 https://github.com/bio15anu/revelio.git revelio_repo

    # Add MD tags using samtools calmd (required for revelio)
    samtools calmd -b ${bam} ${fasta} > ${prefix}.calmd.bam 2>/dev/null

    # Run revelio to mask bisulfite conversions
    python revelio_repo/revelio.py ${prefix}.calmd.bam ${prefix}.revelio.bam

    # Index the output BAM
    samtools index ${prefix}.revelio.bam

    # Clean up intermediate files
    rm -f ${prefix}.calmd.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revelio: \$(cd revelio_repo && git rev-parse --short HEAD)
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
