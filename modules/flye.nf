/*
 * FLYE - Long-read metagenomic assembly
 *
 * Assembles ONT reads in metagenome mode. Replaces SPAdes for
 * long-read AMR subset assembly. Handles empty input gracefully.
 */

process FLYE {

    label 'flye'
    container 'quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0'

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("${sample}_contigs.fasta"), emit: contigs
        path("versions.yml"),                                emit: versions

    script:
    """
    # Check if input reads are non-empty
    READ_COUNT=\$(zcat ${reads} 2>/dev/null | head -1 | wc -l)

    if [ -s ${reads} ] && [ "\$READ_COUNT" -gt 0 ]; then
        flye \\
            ${params.flye_read_type} ${reads} \\
            --meta \\
            --out-dir flye_out \\
            --threads ${task.cpus} || echo "Flye failed to assemble, creating empty fasta"

        if [ -f flye_out/assembly.fasta ]; then
            cp flye_out/assembly.fasta ${sample}_contigs.fasta
        else
            touch ${sample}_contigs.fasta
        fi
    else
        # Empty input - skip assembly
        touch ${sample}_contigs.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
