/*
 * BWA_MEM_SNP - Align paired-end short reads to CDS reference for SNP profiling
 *
 * Aligns reads using BWA-MEM, pipes through samtools for MAPQ filtering,
 * sorting, and indexing. Same piped pattern as spike_in_removal.nf.
 *
 * Output BAM is used by EXTRACT_CODON_FREQS for per-position AA frequency analysis.
 */

process BWA_MEM_SNP {

    label 'lowmem'
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'

    input:
        tuple val(sample), path(r1), path(r2)
        path(bwa_index)    // directory from BWA_INDEX_SNP

    output:
        tuple val(sample), path("${sample}_snp.bam"), path("${sample}_snp.bam.bai"), emit: bam
        path("${sample}_snp_align_stats.txt"), emit: stats
        path("versions.yml"),                  emit: versions

    script:
    def idx_base = "${bwa_index}/reference.fa"
    """
    set -o pipefail

    bwa mem \\
        -t ${task.cpus} \\
        ${idx_base} \\
        ${r1} \\
        ${r2} \\
        2> ${sample}_snp_align_stats.txt \\
    | samtools view -bS -q ${params.snp_min_mapq} - \\
    | samtools sort -@ ${task.cpus} -o ${sample}_snp.bam -

    samtools index ${sample}_snp.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -m1 'Version' | sed 's/Version: //' || echo "unknown")
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    END_VERSIONS
    """
}
