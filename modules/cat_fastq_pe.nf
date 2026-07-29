/*
 * CAT_FASTQ_PE - Concatenate paired-end FASTQ files from re-sequencing runs
 *
 * When a sample is sequenced across multiple runs, this process
 * concatenates the R1 and R2 files. Gzip streams are concatenatable,
 * so simple `cat` produces valid gzipped output.
 */

process CAT_FASTQ_PE {

    tag "$sample"
    label 'ultralow'
    container 'ubuntu:22.04'

    input:
        tuple val(sample), path(r1_files), path(r2_files)

    output:
        tuple val(sample), path("${sample}_R1_cat.fastq.gz"), path("${sample}_R2_cat.fastq.gz"), emit: reads
        path("versions.yml"), emit: versions

    script:
    """
    cat ${r1_files.join(' ')} > ${sample}_R1_cat.fastq.gz
    cat ${r2_files.join(' ')} > ${sample}_R2_cat.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat_fastq_pe: "1.0"
    END_VERSIONS
    """
}
