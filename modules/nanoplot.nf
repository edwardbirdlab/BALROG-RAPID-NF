/*
 * NANOPLOT - QC visualization for ONT long reads
 *
 * Generates read length distributions, quality plots, and summary
 * statistics. Emits NanoStats.txt for MultiQC nanostat module.
 */

process NANOPLOT {

    label 'nanoplot'
    container 'quay.io/biocontainers/nanoplot:1.46.2--pyhdfd78af_1'

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("${sample}_nanoplot"), emit: reports
        path("${sample}_nanoplot/${sample}_NanoStats.txt"), emit: stats
        path("versions.yml"),                          emit: versions

    script:
    """
    NanoPlot \\
        --fastq ${reads} \\
        --outdir ${sample}_nanoplot \\
        --prefix ${sample}_ \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed -e 's/NanoPlot //' || echo "unknown")
    END_VERSIONS
    """
}
