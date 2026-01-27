process MERGE_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.6' :
        'biocontainers/pigz:2.6' }"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        cat ${reads} > ${prefix}.merged.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //')
        END_VERSIONS
        """
    } else {
        // For paired-end, reads come as [ [r1_lane1, r1_lane2], [r2_lane1, r2_lane2] ]
        // or flattened. We need to separate R1 and R2 files
        def r1_files = reads.findAll { it.toString().contains('_R1') || it.toString().contains('_1.') }
        def r2_files = reads.findAll { it.toString().contains('_R2') || it.toString().contains('_2.') }
        """
        cat ${r1_files.join(' ')} > ${prefix}_R1.merged.fastq.gz
        cat ${r2_files.join(' ')} > ${prefix}_R2.merged.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}.merged.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //')
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}_R1.merged.fastq.gz
        touch ${prefix}_R2.merged.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //')
        END_VERSIONS
        """
    }
}
