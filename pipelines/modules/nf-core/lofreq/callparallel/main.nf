process LOFREQ_CALLPARALLEL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4'
        : 'biocontainers/lofreq:2.1.5--py38h588ecb2_4'}"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_intervals = intervals ? "-l ${intervals}" : ""

    def alignment_cram = bam.Extension == "cram" ? true : false
    def alignment_out = alignment_cram ? bam.BaseName + ".bam" : "${bam}"

    def samtools_cram_convert = ''
    samtools_cram_convert += alignment_cram ? "    samtools view -T ${fasta} ${bam} -@ ${task.cpus} -o ${alignment_out}\n" : ''
    samtools_cram_convert += alignment_cram ? "    samtools index ${alignment_out}\n" : ''

    def samtools_cram_remove = ''
    samtools_cram_remove += alignment_cram ? "    rm ${alignment_out}\n" : ''
    samtools_cram_remove += alignment_cram ? "    rm ${alignment_out}.bai\n " : ''
    """
    ${samtools_cram_convert}

    # Check if BAM has sufficient aligned reads for variant calling
    MAPPED_READS=\$(samtools view -c -F 4 ${alignment_out})
    echo "Mapped reads in BAM: \${MAPPED_READS}"

    if [ "\${MAPPED_READS}" -lt 100 ]; then
        echo "WARNING: BAM has insufficient aligned reads (\${MAPPED_READS}) for variant calling"
        echo "Creating empty VCF file to allow pipeline continuation"

        # Create empty VCF with proper header
        echo '##fileformat=VCFv4.2' > ${prefix}.vcf
        echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">' >> ${prefix}.vcf
        echo '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">' >> ${prefix}.vcf
        echo '##FILTER=<ID=PASS,Description="All filters passed">' >> ${prefix}.vcf
        echo '##contig=<ID=chr1>' >> ${prefix}.vcf
        echo -e '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' >> ${prefix}.vcf
        bgzip ${prefix}.vcf
        tabix -p vcf ${prefix}.vcf.gz
    else
        # Run LoFreq variant calling
        # NOTE: LoFreq has HARDCODED /tmp/ paths for Bonferroni correction files
        # that ignore TMPDIR. On Fusion, this causes race conditions.
        # FIX: awsbatch.config mounts NVMe at /tmp via: volumes = '/scratch/fusion:/tmp'
        lofreq \\
            call-parallel \\
            --pp-threads ${task.cpus} \\
            ${args} \\
            ${options_intervals} \\
            -f ${fasta} \\
            -o ${prefix}.vcf.gz \\
            ${alignment_out}
    fi

    ${samtools_cram_remove}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
