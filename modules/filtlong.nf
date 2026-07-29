/*
 * FILTLONG - Quality and length filtering for long reads
 *
 * Filters ONT reads by minimum length and keeps the best reads
 * by quality score up to the specified percentage.
 */

process FILTLONG {

    label 'lowmem'
    container 'quay.io/biocontainers/filtlong:0.3.1--h077b44d_0'

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("${sample}_filtered.fastq.gz"), emit: filtered
        path("${sample}_filtlong.log"),                          emit: log
        path("versions.yml"),                                    emit: versions

    script:
    """
    filtlong \\
        --min_length ${params.filtlong_min_length} \\
        --keep_percent ${params.filtlong_keep_percent} \\
        ${reads} 2> ${sample}_filtlong.log | gzip > ${sample}_filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filtlong: \$(filtlong --version 2>&1 | sed -e 's/Filtlong v//' || echo "unknown")
    END_VERSIONS
    """
}
