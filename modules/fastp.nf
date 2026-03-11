process FASTP {

    label 'lowmem'
    container 'biocontainers/fastp:v0.20.1_cv1'

    input:
        tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), path("${sample}_trimmed_R1.fastq.gz"), path("${sample}_trimmed_R2.fastq.gz"), emit: trimmed_fastq
        path("${sample}_fastp.json"), emit: json
        path("${sample}_fastp.html"), emit: html
        path("versions.yml"),         emit: versions

    script:
    """
    fastp \
        -i ${r1} -I ${r2} \
        -o ${sample}_trimmed_R1.fastq.gz \
        -O ${sample}_trimmed_R2.fastq.gz \
        -q ${params.fastp_q} \
        -l ${params.fastp_minlen} \
        --thread ${task.cpus} \
        -j ${sample}_fastp.json \
        -h ${sample}_fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //' || echo "unknown")
    END_VERSIONS
    """
}
