/*
 * EXTRACT_BACTERIAL_READS_SE - KrakenTools bacterial read extraction (single-end)
 *
 * Extracts reads classified as Bacteria (taxid 2, including all children)
 * from Kraken2 output. Used to pre-filter reads before Nonpareil coverage
 * estimation, removing host contamination signal that would bias results.
 *
 * Generates single-end gzipped FASTQ and a stats TSV summarizing extraction.
 */

process EXTRACT_BACTERIAL_READS_SE {

    label 'ultralow'
    container 'quay.io/biocontainers/krakentools:1.2.1--pyh7e72e81_0'

    input:
        tuple val(sample), path(reads), path(k2_output), path(k2_report)

    output:
        tuple val(sample), path("${sample}_bacteria.fastq.gz"), emit: bacterial_reads
        path("${sample}_bacterial_extraction_stats.tsv"),        emit: stats
        path("versions.yml"),                                    emit: versions

    script:
    """
    # Extract reads classified under Bacteria (taxid 2) including all descendant taxa
    extract_kraken_reads.py \\
        -k ${k2_output} \\
        -r ${k2_report} \\
        -s ${reads} \\
        -o ${sample}_bacteria.fastq \\
        -t 2 \\
        --include-children \\
        --fastq-output

    # Compress output (extract_kraken_reads.py writes uncompressed FASTQ)
    gzip ${sample}_bacteria.fastq

    # Generate extraction statistics
    TOTAL=\$(zcat ${reads} | awk 'END{print NR/4}')
    BACT=\$(zcat ${sample}_bacteria.fastq.gz | awk 'END{print NR/4}')

    if [ "\${TOTAL}" -gt 0 ]; then
        PCT=\$(awk "BEGIN {printf \\"%.2f\\", (\${BACT}/\${TOTAL})*100}")
    else
        PCT="0.00"
    fi

    echo -e "sample\\ttotal_reads\\tbacterial_reads\\tbacterial_pct" > ${sample}_bacterial_extraction_stats.tsv
    echo -e "${sample}\\t\${TOTAL}\\t\${BACT}\\t\${PCT}" >> ${sample}_bacterial_extraction_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: \$(echo "1.2.1")
    END_VERSIONS
    """
}
