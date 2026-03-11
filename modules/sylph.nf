process SYLPH {

    label 'sylph_profile'
    container 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'

    input:
        tuple val(sample), path(r1), path(r2)
        path(db)

    output:
        tuple val(sample), path("${sample}_profile.tsv"), emit: profile
        path("versions.yml"),                              emit: versions

    script:
    """
    sylph profile ${db} \
        -o ${sample}_profile.tsv \
        -t ${task.cpus} \
        ${r1} ${r2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
