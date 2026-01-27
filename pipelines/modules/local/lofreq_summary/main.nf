process LOFREQ_SUMMARY {
    tag "LoFreq variant summary"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/astral-sh/uv:python3.12-bookworm':
        'ghcr.io/astral-sh/uv:python3.12-bookworm' }"

    input:
    path vcf_files
    path(script)  // Staged script for cloud execution

    output:
    path "lofreq_mqc.json", emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Create temp directory for VCF files
    mkdir -p vcf_temp

    # Copy VCF files to temp directory to organize them
    for vcf in ${vcf_files}; do
        cp "\$vcf" vcf_temp/
    done

    # Run the summary script
    uv run ${script} \\
        --input-dir vcf_temp \\
        --output lofreq_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo '{"id": "lofreq_variants", "section_name": "LoFreq Variant Calling", "data": {}}' > lofreq_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}