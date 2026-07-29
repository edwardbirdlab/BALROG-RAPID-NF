/*
 * BWA_INDEX_SNP - Index CDS FASTA for targeted SNP profiling
 *
 * Builds a BWA index from the user-provided CDS nucleotide FASTA.
 * Runs once per pipeline execution (not per sample).
 * The index is reused across all short-read samples.
 */

process BWA_INDEX_SNP {

    label 'ultralow'
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'

    input:
        path(cds_fasta)

    output:
        path("bwa_index"), emit: index
        path("versions.yml"), emit: versions

    script:
    """
    mkdir bwa_index
    cp ${cds_fasta} bwa_index/reference.fa
    bwa index bwa_index/reference.fa
    samtools faidx bwa_index/reference.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -m1 'Version' | sed 's/Version: //' || echo "unknown")
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo "unknown")
    END_VERSIONS
    """
}
