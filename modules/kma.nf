/*
 * KMA - Align reads to a COI reference database for insect ID
 *
 * Emits .res (per-template coverage/identity/depth) for SUMMARIZE_COI.
 * -mem_mode and -and are load-bearing: without them KMA produced no .res
 * at all in testing. Exit code 95 is a documented non-error KMA condition
 * (see nf-core/modules#7251); treat it as success.
 *
 * The index prefix is auto-detected from whatever *.comp.b file is in the
 * given directory, rather than assumed -- different KMA databases use
 * different prefix conventions (e.g. bold_coi_diptera_v5.comp.b).
 */

process KMA {

    label 'kma'
    // Digest-pinned, no tag: Apptainer/Singularity rejects a combined
    // tag@digest reference ("Docker references with both a tag and digest
    // are currently not supported"), even though Docker accepts it. The
    // digest alone is still fully pinned -- this is kma:1.6.13--h118bc1c_0.
    container 'quay.io/biocontainers/kma@sha256:90e62ef87fd8ff3bdcb2ce316dfa69faf854396e7e3a18bd57d7a2b08920759d'

    input:
        tuple val(sample), path(r1), path(r2)
        path(kma_db)   // directory containing the KMA index (<prefix>.comp.b, .length.b, .seq.b, .name)

    output:
        tuple val(sample), path("${sample}.res"),     emit: res
        tuple val(sample), path("${sample}.mapstat"), emit: mapstat
        path("versions.yml"),                          emit: versions

    script:
    """
    INDEX_FILE=\$(ls ${kma_db}/*.comp.b | head -n 1)
    if [ -z "\$INDEX_FILE" ]; then
        echo "ERROR: no *.comp.b index file found in ${kma_db}" >&2
        exit 1
    fi
    INDEX_BASE=\${INDEX_FILE%.comp.b}

    kma -ipe ${r1} ${r2} \\
        -o ${sample} \\
        -t_db \$INDEX_BASE \\
        -t ${task.cpus} \\
        -1t1 -mem_mode -and -apm f -ef \\
        || [ \$? -eq 95 ]

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(kma -v 2>&1 | sed 's/KMA-//' || echo "unknown")
    END_VERSIONS
    """
}
