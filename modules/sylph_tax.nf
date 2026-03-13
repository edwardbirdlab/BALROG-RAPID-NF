/*
 * SYLPH_TAX
 *
 * Converts raw sylph profile TSV output into taxonomic profiles (.sylphmpa)
 * using sylph-tax taxprof. The .sylphmpa files are natively parsed by
 * MultiQC's sylphtax module for taxonomy bargraphs and general stats.
 *
 * Requires a pre-downloaded sylph-tax taxonomy database directory
 * (created via: sylph-tax download --download-to /path/to/sylph_tax_db)
 * with a config.json inside pointing taxonomy_dir to "."
 */

process SYLPH_TAX {

    label 'ultralow'
    container 'quay.io/biocontainers/sylph-tax:1.8.0--pyhdfd78af_0'

    input:
        tuple val(sample), path(profile_tsv)
        path(taxonomy_db)
        val(taxonomy_name)

    output:
        path("*.sylphmpa"), emit: sylphmpa
        path("versions.yml"), emit: versions

    script:
    """
    # Point sylph-tax at the staged taxonomy database directory
    export SYLPH_TAXONOMY_CONFIG="\${PWD}/${taxonomy_db}/config.json"

    sylph-tax taxprof ${profile_tsv} -t ${taxonomy_name}

    # If no .sylphmpa was produced (e.g., no hits), create an empty one
    if ! ls *.sylphmpa 1>/dev/null 2>&1; then
        echo -e "clade_name\\trelative_abundance\\tsequence_abundance\\tANI" > ${sample}_empty.sylphmpa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
