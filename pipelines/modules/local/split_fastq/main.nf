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
            -j ${task.cpus}

        # Rename and gzip outputs (seqkit outputs as xxx.read1.fastq, xxx.read2.fastq)
        for f in split_output/*.read1.fastq; do
            part=\$(basename \$f .read1.fastq)
            gzip -c \$f > ${prefix}_R1.part_\${part}.fastq.gz
        done
        for f in split_output/*.read2.fastq; do
            part=\$(basename \$f .read2.fastq)
            gzip -c \$f > ${prefix}_R2.part_\${part}.fastq.gz
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
            -j ${task.cpus}

        # Rename and gzip outputs (seqkit outputs as 001.fastq, 002.fastq, etc.)
        for f in split_output/*.fastq; do
            part=\$(basename \$f .fastq)
            gzip -c \$f > ${prefix}_R1.part_\${part}.fastq.gz
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | sed 's/seqkit v//')
        END_VERSIONS
        """
    }
}
