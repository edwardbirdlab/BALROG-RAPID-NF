process SPIKE_IN_REMOVAL {

    label 'lowmem'
    container 'biocontainers/bowtie2:v2.4.1_cv1'

    input:
        tuple val(sample), path(r1), path(r2)
        path(bt2_index)    // directory containing *.bt2 files

    output:
        tuple val(sample), path("${sample}_spike_removed_R1.fastq.gz"), path("${sample}_spike_removed_R2.fastq.gz"), emit: cleaned_reads
        path("${sample}_bowtie2_spike.log"), emit: log
        path("${sample}_spike_stats.tsv"),   emit: stats
        path("versions.yml"),                emit: versions

    script:
    def idx_base = "${bt2_index}/t_thermophilus"
    """
    bowtie2 \\
        -p ${task.cpus} \\
        -x ${idx_base} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --very-sensitive-local \\
        --un-conc-gz ${sample}_spike_removed \\
        > /dev/null \\
        2> ${sample}_bowtie2_spike.log

    # Bowtie2 --un-conc-gz creates {prefix}.1 and {prefix}.2
    mv ${sample}_spike_removed.1 ${sample}_spike_removed_R1.fastq.gz
    mv ${sample}_spike_removed.2 ${sample}_spike_removed_R2.fastq.gz

    # Parse alignment stats for easy downstream use
    TOTAL=\$(grep "reads; of these:" ${sample}_bowtie2_spike.log | awk '{print \$1}')
    ALIGNED=\$(grep "aligned concordantly" ${sample}_bowtie2_spike.log | head -1 | awk '{print \$1}')
    RATE=\$(grep "overall alignment rate" ${sample}_bowtie2_spike.log | awk '{print \$1}')
    echo -e "sample\\ttotal_reads\\tspike_in_reads\\talignment_rate" > ${sample}_spike_stats.tsv
    echo -e "${sample}\\t\${TOTAL}\\t\${ALIGNED}\\t\${RATE}" >> ${sample}_spike_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/.*version //' || echo "unknown")
    END_VERSIONS
    """
}
