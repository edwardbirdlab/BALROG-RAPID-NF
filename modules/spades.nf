process SPADES {

    label 'spades_subset'
    container 'ebird013/spades:3.15.5'

    input:
        tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), path("${sample}_contigs.fasta"), emit: contigs
        path("versions.yml"),                                emit: versions

    script:
    """
    # Check if input reads are non-empty (have at least one read)
    READ_COUNT=\$(zcat ${r1} 2>/dev/null | head -1 | wc -l)

    if [ -s ${r1} ] && [ "\$READ_COUNT" -gt 0 ]; then
        spades.py \
            -1 ${r1} \
            -2 ${r2} \
            -o spades_out \
            --threads ${task.cpus} \
            --memory ${task.memory.toGiga()} \
            --only-assembler \
            --meta

        if [ -f spades_out/contigs.fasta ]; then
            cp spades_out/contigs.fasta ${sample}_contigs.fasta
        else
            touch ${sample}_contigs.fasta
        fi
    else
        # Empty input - skip assembly
        touch ${sample}_contigs.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
