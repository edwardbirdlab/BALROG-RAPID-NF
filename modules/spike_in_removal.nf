/*
 * SPIKE_IN_REMOVAL - Bowtie2 alignment + samtools extraction (piped)
 *
 * Aligns reads to T. thermophilus spike-in reference, then extracts
 * unmapped read pairs (both mates unmapped) directly via pipe.
 * No intermediate SAM/BAM file is written to disk.
 *
 * Bowtie2 stats are captured from stderr for MultiQC integration.
 */

process SPIKE_IN_REMOVAL {

    label 'lowmem'
    container 'quay.io/biocontainers/mulled-v2-229691629e0b12c862d76101f90a597d5c1c81d4:484c804e1d5952c9023891b6f9a19f7f15815145-0'

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
    set -o pipefail

    # Pipe: bowtie2 → samtools view (unmapped pairs) → sort → fastq
    # bowtie2 stderr = alignment stats (for MultiQC), stdout = SAM records
    # -f 12:  both read and mate unmapped
    # -F 256: exclude secondary alignments
    bowtie2 \\
        -p ${task.cpus} \\
        -x ${idx_base} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --very-sensitive \\
        2> ${sample}_bowtie2_spike.log \\
    | samtools view -b -f 12 -F 256 - \\
    | samtools sort -n -m ${task.memory.toGiga()}G -@ ${task.cpus} - \\
    | samtools fastq -@ ${task.cpus} - \\
        -1 ${sample}_spike_removed_R1.fastq.gz \\
        -2 ${sample}_spike_removed_R2.fastq.gz

    # Parse spike-in stats from bowtie2 log (stderr)
    TOTAL=\$(grep -m1 'reads; of these' ${sample}_bowtie2_spike.log | awk '{print \$1}' || true)
    RATE=\$(grep 'overall alignment rate' ${sample}_bowtie2_spike.log | awk '{print \$1}' || true)

    # Guard against empty input (0 reads)
    if [ -z "\${TOTAL}" ] || [ "\${TOTAL}" -eq 0 ] 2>/dev/null; then
        TOTAL=0
        MAPPED=0
        RATE="0.00%"
    else
        RATE_NUM=\$(echo "\${RATE}" | sed 's/%//')
        MAPPED=\$(awk "BEGIN {printf \\"%d\\", (\${RATE_NUM}/100)*\${TOTAL}}")
    fi

    echo -e "sample\\ttotal_read_pairs\\tspike_in_pairs\\talignment_rate" > ${sample}_spike_stats.tsv
    echo -e "${sample}\\t\${TOTAL}\\t\${MAPPED}\\t\${RATE}" >> ${sample}_spike_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/.*version //' || echo "unknown")
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    END_VERSIONS
    """
}
