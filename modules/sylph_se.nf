/*
 * SYLPH_SE - Single-end Sylph profiling
 *
 * For single-end reads, sylph requires an explicit sketch step
 * before profiling, unlike paired-end which auto-sketches.
 */

process SYLPH_SE {

    label 'sylph_profile'
    container 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'

    input:
        tuple val(sample), path(reads)
        path(db)

    output:
        tuple val(sample), path("${sample}_profile.tsv"), emit: profile
        path("versions.yml"),                              emit: versions

    script:
    """
    
    sylph profile ${db} \
        -o ${sample}_profile.tsv \
        -t ${task.cpus} \
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
