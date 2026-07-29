/*
 * CAT_FASTQ_SE - Concatenate single-end FASTQ files from re-sequencing runs
 *
 * When a long-read sample is sequenced across multiple runs, this process
 * concatenates the read files. Gzip streams are concatenatable,
 * so simple `cat` produces valid gzipped output.
 */

process CAT_FASTQ_SE {

    tag "$sample"
    label 'ultralow'
    container 'ubuntu:22.04'

    input:
        tuple val(sample), path(read_files)

    output:
        tuple val(sample), path("${sample}_cat.fastq.gz"), emit: reads
        path("versions.yml"), emit: versions

    script:
    """
    cat ${read_files.join(' ')} > ${sample}_cat.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat_fastq_se: "1.0"
    END_VERSIONS
    """
}
