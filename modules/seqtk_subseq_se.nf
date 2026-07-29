process SEQTK_SUBSEQ_SE {

    label 'ultralow'
    container 'quay.io/biocontainers/seqtk:1.5--h577a1d6_1'

    input:
        tuple val(sample), path(reads), path(read_ids)

    output:
        tuple val(sample), path("${sample}_amr_reads.fastq.gz"), emit: subset_reads
        path("versions.yml"),                                     emit: versions

    script:
    """
    if [ -s ${read_ids} ]; then
        seqtk subseq ${reads} ${read_ids} | gzip > ${sample}_amr_reads.fastq.gz
    else
        # No Diamond hits - create empty gzipped file
        echo -n | gzip > ${sample}_amr_reads.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | head -3 | tail -1 || echo "unknown")
    END_VERSIONS
    """
}
