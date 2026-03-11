process SEQTK_SUBSEQ {

    label 'ultralow'
    container 'quay.io/biocontainers/seqtk:1.5--h577a1d6_1'

    input:
        tuple val(sample), path(r1), path(r2), path(read_ids)

    output:
        tuple val(sample), path("${sample}_amr_R1.fastq.gz"), path("${sample}_amr_R2.fastq.gz"), emit: subset_reads
        path("versions.yml"),                                                                      emit: versions

    script:
    """
    if [ -s ${read_ids} ]; then
        seqtk subseq ${r1} ${read_ids} | gzip > ${sample}_amr_R1.fastq.gz
        seqtk subseq ${r2} ${read_ids} | gzip > ${sample}_amr_R2.fastq.gz
    else
        # No Diamond hits - create empty gzipped files
        echo -n | gzip > ${sample}_amr_R1.fastq.gz
        echo -n | gzip > ${sample}_amr_R2.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | head -3 | tail -1 || echo "unknown")
    END_VERSIONS
    """
}
