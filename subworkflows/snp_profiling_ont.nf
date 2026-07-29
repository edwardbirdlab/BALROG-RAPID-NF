/*
 * SNP Profiling Subworkflow (Long-read, Single-end, ONT)
 *
 * Targeted amino acid variant profiling at user-specified CDS positions:
 *   1. minimap2 alignment of SE long reads to CDS reference (per sample)
 *   2. Extract per-read codon frequencies at target positions (per sample)
 *   3. Summarize for MultiQC custom content (per sample)
 *
 * Uses minimap2 instead of BWA-MEM for long reads. minimap2 indexes inline
 * (fast for small CDS references), so no separate index step is needed.
 *
 * Designed for pooled/metagenomic samples: reports AA frequency distributions
 * rather than consensus genotypes.
 */

include { MINIMAP2_SNP                                       } from '../modules/minimap2_snp'
include { EXTRACT_CODON_FREQS   as EXTRACT_CODON_FREQS_SE    } from '../modules/extract_codon_freqs'
include { SUMMARIZE_SNP_PROFILING as SUMMARIZE_SNP_PROFILING_ONT } from '../modules/summarize_snp_profiling'


workflow SNP_PROFILING_ONT {

    take:
        ch_reads          // tuple(sample, reads)
        ch_cds_fasta      // path to CDS nucleotide FASTA
        ch_positions_csv  // path to positions of interest CSV

    main:
        // Step 1: Align long reads to CDS reference (minimap2 indexes inline)
        MINIMAP2_SNP(ch_reads, ch_cds_fasta)

        // Step 2: Extract codon frequencies at target positions
        EXTRACT_CODON_FREQS_SE(MINIMAP2_SNP.out.bam, ch_cds_fasta, ch_positions_csv)

        // Step 3: Generate MultiQC custom content files
        SUMMARIZE_SNP_PROFILING_ONT(EXTRACT_CODON_FREQS_SE.out.codon_freqs)

    emit:
        multiqc_snp = SUMMARIZE_SNP_PROFILING_ONT.out.generalstats
            .mix(SUMMARIZE_SNP_PROFILING_ONT.out.detail)
        codon_freqs = EXTRACT_CODON_FREQS_SE.out.codon_freqs
        versions    = MINIMAP2_SNP.out.versions
            .mix(EXTRACT_CODON_FREQS_SE.out.versions,
                 SUMMARIZE_SNP_PROFILING_ONT.out.versions)
}
