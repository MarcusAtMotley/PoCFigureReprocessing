process SPLIT_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(chunk_size)

    output:
    tuple val(meta), path("*.part_*.bam"), path("*.part_*.bam.bai"), emit: bams
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    # Get total read count
    total_reads=\$(samtools view -c ${bam})

    # Calculate number of chunks
    num_chunks=\$(( (total_reads + ${chunk_size} - 1) / ${chunk_size} ))

    # If only 1 chunk needed, just copy the file
    if [ \$num_chunks -le 1 ]; then
        cp ${bam} ${prefix}.part_001.bam
        cp ${bai} ${prefix}.part_001.bam.bai
    else
        # Extract header
        samtools view -H ${bam} > header.sam

        # Split reads into chunks
        samtools view ${bam} | split -l ${chunk_size} -d -a 3 --additional-suffix=.sam - chunk_

        # Convert each chunk to BAM with header
        for chunk in chunk_*.sam; do
            part=\$(echo \$chunk | sed 's/chunk_//' | sed 's/.sam//')
            # Pad part number
            padded=\$(printf "%03d" \$((10#\$part + 1)))

            cat header.sam \$chunk | samtools view -b -o ${prefix}.part_\${padded}.bam -
            samtools index ${prefix}.part_\${padded}.bam
            rm \$chunk
        done

        rm header.sam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
