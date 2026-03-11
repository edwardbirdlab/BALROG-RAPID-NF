process DIAMOND_BLASTX {

    label 'diamond'
    container 'quay.io/biocontainers/diamond:2.1.24--hf93d47f_0'

    input:
        tuple val(sample), path(r1), path(r2)
        path(db)

    output:
        tuple val(sample), path("${sample}_diamond_hits.tsv"),  emit: hits
        tuple val(sample), path("${sample}_hit_readids.txt"),   emit: read_ids
        path("versions.yml"),                                    emit: versions

    script:
    """
    diamond blastx \
        --db ${db} \
        --query ${r1} \
        --out ${sample}_diamond_hits.tsv \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue ${params.diamond_evalue} \
        --max-target-seqs ${params.diamond_max_targets} \
        --threads ${task.cpus}

    # Extract unique read IDs for seqtk subseq
    if [ -s ${sample}_diamond_hits.tsv ]; then
        cut -f1 ${sample}_diamond_hits.tsv | sort -u > ${sample}_hit_readids.txt
    else
        touch ${sample}_hit_readids.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
