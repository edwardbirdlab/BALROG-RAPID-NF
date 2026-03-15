/*
 * EXTRACT_BACTERIAL_READS - KrakenTools bacterial read extraction
 *
 * Extracts reads classified as Bacteria (taxid 2, including all children)
 * from Kraken2 output. Used to pre-filter reads before Nonpareil coverage
 * estimation, removing host contamination signal that would bias results.
 *
 * Generates paired-end gzipped FASTQ and a stats TSV summarizing extraction.
 */

process EXTRACT_BACTERIAL_READS {

    label 'ultralow'
    container 'quay.io/biocontainers/krakentools:1.2.1--pyh7e72e81_0'

    input:
        tuple val(sample), path(r1), path(r2), path(k2_output), path(k2_report)

    output:
        tuple val(sample), path("${sample}_bacteria_R1.fastq.gz"), path("${sample}_bacteria_R2.fastq.gz"), emit: bacterial_reads
        path("${sample}_bacterial_extraction_stats.tsv"), emit: stats
        path("versions.yml"),                             emit: versions

    script:
    """
    # Extract reads classified under Bacteria (taxid 2) including all descendant taxa
    extract_kraken_reads.py \\
        -k ${k2_output} \\
        -r ${k2_report} \\
        -s ${r1} \\
        -s2 ${r2} \\
        -o ${sample}_bacteria_R1.fastq \\
        -o2 ${sample}_bacteria_R2.fastq \\
        -t 2 \\
        --include-children \\
        --fastq-output

    # Compress outputs (extract_kraken_reads.py writes uncompressed FASTQ)
    gzip ${sample}_bacteria_R1.fastq
    gzip ${sample}_bacteria_R2.fastq

    # Generate extraction statistics
    TOTAL_R1=\$(zcat ${r1} | awk 'END{print NR/4}')
    BACT_R1=\$(zcat ${sample}_bacteria_R1.fastq.gz | awk 'END{print NR/4}')

    if [ "\${TOTAL_R1}" -gt 0 ]; then
        PCT=\$(awk "BEGIN {printf \\"%.2f\\", (\${BACT_R1}/\${TOTAL_R1})*100}")
    else
        PCT="0.00"
    fi

    echo -e "sample\\ttotal_read_pairs\\tbacterial_read_pairs\\tbacterial_pct" > ${sample}_bacterial_extraction_stats.tsv
    echo -e "${sample}\\t\${TOTAL_R1}\\t\${BACT_R1}\\t\${PCT}" >> ${sample}_bacterial_extraction_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: \$(echo "1.2.1")
    END_VERSIONS
    """
}
