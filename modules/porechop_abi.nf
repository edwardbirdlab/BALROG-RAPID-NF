/*
 * PORECHOP_ABI - ONT adapter trimming with ab initio adapter discovery
 *
 * Porechop_ABI discovers adapters de novo from the reads rather than
 * relying on a predefined adapter list. This makes it more robust for
 * novel or unexpected adapter contamination in ONT libraries.
 */

process PORECHOP_ABI {

    label 'porechop'
    container 'quay.io/biocontainers/porechop_abi:0.5.1--py311h2de2dd3_0'

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("${sample}_pc.fastq.gz"), emit: trimmed
        path("${sample}_porechop.log"),                    emit: log
        path("versions.yml"),                              emit: versions

    script:
    """
    porechop_abi \\
        -abi \\
        -i ${reads} \\
        -o ${sample}_pc.fastq.gz \\
        --threads ${task.cpus} \\
        2>&1 | tee ${sample}_porechop.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop_abi: \$(porechop_abi --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
