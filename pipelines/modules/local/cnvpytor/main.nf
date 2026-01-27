process CNVPYTOR_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cnvpytor=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.3.1--pyhdfd78af_0' :
        'biocontainers/cnvpytor:1.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.pytor")    , emit: pytor
    tuple val(meta), path("*_cnv.tsv")  , emit: calls
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bin_size = params.cnv_bin_size ?: 10000
    """
    # Create CNVpytor file from BAM
    cnvpytor -root ${prefix}.pytor -rd ${bam}

    # Calculate histogram with specified bin size
    cnvpytor -root ${prefix}.pytor -his ${bin_size}

    # Partition the data
    cnvpytor -root ${prefix}.pytor -partition ${bin_size}

    # Call CNVs
    cnvpytor -root ${prefix}.pytor -call ${bin_size} > ${prefix}_cnv.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version 2>&1 | sed 's/CNVpytor //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pytor
    touch ${prefix}_cnv.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version 2>&1 | sed 's/CNVpytor //')
    END_VERSIONS
    """
}
