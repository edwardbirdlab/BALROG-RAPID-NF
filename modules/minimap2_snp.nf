/*
 * MINIMAP2_SNP - Align single-end long reads to CDS reference for SNP profiling
 *
 * Aligns ONT reads using minimap2 (indexes inline — fast for small CDS references),
 * pipes through samtools for MAPQ filtering, sorting, and indexing.
 *
 * Output BAM is used by EXTRACT_CODON_FREQS for per-position AA frequency analysis.
 */

process MINIMAP2_SNP {

    label 'lowmem'
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:e1ea28074233d7265a5dc2111d6e55130dff5653-2'

    input:
        tuple val(sample), path(reads)
        path(cds_fasta)

    output:
        tuple val(sample), path("${sample}_snp.bam"), path("${sample}_snp.bam.bai"), emit: bam
        path("${sample}_snp_align_stats.txt"), emit: stats
        path("versions.yml"),                  emit: versions

    script:
    """
    set -o pipefail

    minimap2 \\
        -a \\
        -x map-ont \\
        -t ${task.cpus} \\
        ${cds_fasta} \\
        ${reads} \\
        2> ${sample}_snp_align_stats.txt \\
    | samtools view -bS -q ${params.snp_min_mapq} - \\
    | samtools sort -@ ${task.cpus} -o ${sample}_snp.bam -

    samtools index ${sample}_snp.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1 || echo "unknown")
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    END_VERSIONS
    """
}
