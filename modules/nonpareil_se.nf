/*
 * NONPAREIL_SE - Metagenomic coverage estimation (single-end)
 *
 * Estimates sequencing coverage by analyzing read redundancy.
 * Generates JSON output for MultiQC integration.
 */

process NONPAREIL_SE {

    label 'lowmem'
    container 'quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2'

    input:
        tuple val(sample), path(reads)

    output:
        path("${sample}.npo"),              emit: npo
        path("${sample}.npa"),              emit: npa
        path("${sample}_nonpareil.json"),   emit: json
        path("versions.yml"),               emit: versions

    script:
    """
    nonpareil \\
        -s ${reads} \\
        -T kmer \\
        -f fastq \\
        -b ${sample} \\
        -t ${task.cpus}

    # Generate JSON for MultiQC nonpareil module
    NonpareilCurves.R \\
        --json=${sample}_nonpareil.json \\
        --labels=${sample} \\
        ${sample}.npo

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(nonpareil -V 2>&1 | sed 's/Nonpareil v//' || echo "unknown")
    END_VERSIONS
    """
}
