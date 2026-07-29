process FASTQC_SE {

    label 'ultralow'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("${sample}_fastqc"), emit: reports
        path("${sample}_fastqc/*.zip"),              emit: zip
        path("versions.yml"),                        emit: versions

    script:
    """
    mkdir ${sample}_fastqc
    ln -s ${reads} ${sample}.fastq.gz
    fastqc -o ${sample}_fastqc -t ${task.cpus} --memory ${Math.min(task.memory.toMega(), 10000)} ${sample}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version 2>&1 | sed -e 's/FastQC v//' || echo "unknown")
    END_VERSIONS
    """
}
