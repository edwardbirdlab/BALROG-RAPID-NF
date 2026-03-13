/*
 * BBDuk - Adapter/Contaminant Trimming (BBTools Suite)
 *
 * Optional kmer-based trimming for Element Biosciences Aviti reads.
 * Runs after FASTP to remove platform-specific adapter sequences.
 * Enabled with --run_bbduk flag.
 */

process BBDUK {

    label 'lowmem'
    container 'ebird013/bbmap:latest'

    input:
        tuple val(sample), path(r1), path(r2)
        path(adapter_fasta)

    output:
        tuple val(sample), path("${sample}_bbduk_R1.fastq.gz"), path("${sample}_bbduk_R2.fastq.gz"), emit: trimmed
        path("${sample}_bbduk.log"), emit: log
        path("versions.yml"), emit: versions

    script:
    """
    bbduk.sh \\
        in=${r1} \\
        in2=${r2} \\
        out=${sample}_bbduk_R1.fastq.gz \\
        out2=${sample}_bbduk_R2.fastq.gz \\
        ref=${adapter_fasta} \\
        ktrim=${params.bbduk_ktrim} \\
        hdist=${params.bbduk_hdist} \\
        threads=${task.cpus} \\
        ${params.bbduk_additional_args} \\
        2>&1 | tee ${sample}_bbduk.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbmap.sh version 2>&1 | grep "BBMap version" | sed -e "s/BBMap version //g" || echo "unknown")
    END_VERSIONS
    """
}
