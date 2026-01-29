process SPLIT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.5.1--h9ee0642_0' :
        'biocontainers/seqkit:2.5.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    val(chunk_size)

    output:
    tuple val(meta), path("*.part_*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    def is_pe = reads.size() > 1 && !meta.single_end

    if (is_pe) {
        def r1 = reads[0]
        def r2 = reads[1]
        """
        # Split R1 and R2 in sync using seqkit split2
        # When input is gzipped, seqkit outputs gzipped files preserving original name pattern
        seqkit split2 \\
            -1 ${r1} \\
            -2 ${r2} \\
            -s ${chunk_size} \\
            -O split_output \\
            -j ${task.cpus}

        # Move outputs to working directory (seqkit already outputs as .part_NNN.fastq.gz)
        mv split_output/*.part_*.fastq.gz .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | sed 's/seqkit v//')
        END_VERSIONS
        """
    } else {
        def r1 = reads instanceof List ? reads[0] : reads
        """
        # Split SE reads
        seqkit split2 \\
            -1 ${r1} \\
            -s ${chunk_size} \\
            -O split_output \\
            -j ${task.cpus}

        # Move outputs to working directory
        mv split_output/*.part_*.fastq.gz .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | sed 's/seqkit v//')
        END_VERSIONS
        """
    }
}
