process FASTQC {

    label 'ultralow'
    container 'ebird013/fastqc:0.12.1_custom'

    input:
        tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), path("${sample}_fastqc"), emit: reports
        path("${sample}_fastqc/*.zip"),              emit: zip
        path("versions.yml"),                        emit: versions

    script:
    """
    mkdir ${sample}_fastqc
    fastqc -o ${sample}_fastqc -t ${task.cpus} ${r1} ${r2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version 2>&1 | sed -e 's/FastQC v//' || echo "unknown")
    END_VERSIONS
    """
}
