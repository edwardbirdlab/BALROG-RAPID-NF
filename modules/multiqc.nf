process MULTIQC {

    label 'ultralow'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    input:
        path('fastqc_raw/*')
        path('fastp/*')
        path('bbduk/*')
        path('fastqc_trimmed/*')
        path('bowtie2_spike/*')
        path('kraken2_taxonomy/*')
        path('kraken2_host/*')
        path('sylph/*')
        path('amrfinder_summary/*')
        path('nonpareil/*')
        path('custom_qc/*')
        path(multiqc_config)
        path(software_versions)

    output:
        path("multiqc_report.html"), emit: html
        path("multiqc_report_data"), emit: data
        path("versions.yml"),        emit: versions

    script:
    """
    multiqc \\
        -f \\
        -c ${multiqc_config} \\
        --filename multiqc_report.html \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed -e 's/multiqc, version //' || echo "unknown")
    END_VERSIONS
    """
}
