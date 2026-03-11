process KRAKEN2 {

    label 'kraken2'
    container 'ebird013/kraken2:2.1.3'

    input:
        tuple val(sample), path(r1), path(r2)
        path(db)
        val(db_name)

    output:
        tuple val(sample), val(db_name), path("${prefix}_k2report.tsv"), emit: report
        tuple val(sample), val(db_name), path("${prefix}_k2out.tsv"),    emit: output
        path("versions.yml"),                                             emit: versions

    script:
    prefix = db_name ? "${sample}_${db_name}" : "${sample}"
    """
    kraken2 \
        --db ${db} \
        --threads ${task.cpus} \
        --output ${prefix}_k2out.tsv \
        --report ${prefix}_k2report.tsv \
        --minimum-hit-groups ${params.k2_min_hit_groups} \
        --confidence ${params.k2_confidence} \
        --paired ${r1} ${r2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version 2>&1 | head -1 | sed -e 's/Kraken version //' || echo "unknown")
    END_VERSIONS
    """
}
