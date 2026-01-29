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
        seqkit split2 \\
            -1 ${r1} \\
            -2 ${r2} \\
            -s ${chunk_size} \\
            -O split_output \\
            -j ${task.cpus} \\
            --by-size-prefix ""

        # Rename outputs to have consistent naming
        for f in split_output/*_1.part_*.fastq.gz; do
            part=\$(basename \$f | sed 's/.*\\.part_/part_/' | sed 's/_1\\.fastq\\.gz//')
            mv \$f ${prefix}_R1.\${part}.fastq.gz
        done
        for f in split_output/*_2.part_*.fastq.gz; do
            part=\$(basename \$f | sed 's/.*\\.part_/part_/' | sed 's/_2\\.fastq\\.gz//')
            mv \$f ${prefix}_R2.\${part}.fastq.gz
        done

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
            -j ${task.cpus} \\
            --by-size-prefix ""

        # Rename outputs
        for f in split_output/*.part_*.fastq.gz; do
            part=\$(basename \$f | sed 's/.*\\.part_/part_/' | sed 's/\\.fastq\\.gz//')
            mv \$f ${prefix}_R1.\${part}.fastq.gz
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | sed 's/seqkit v//')
        END_VERSIONS
        """
    }
}
