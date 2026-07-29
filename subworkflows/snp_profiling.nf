/*
 * SNP Profiling Subworkflow (Short-read, Paired-end)
 *
 * Targeted amino acid variant profiling at user-specified CDS positions:
 *   1. BWA index of CDS FASTA (once)
 *   2. BWA-MEM alignment of PE reads to CDS reference (per sample)
 *   3. Extract per-read codon frequencies at target positions (per sample)
 *   4. Summarize for MultiQC custom content (per sample)
 *
 * Designed for pooled/metagenomic samples: reports AA frequency distributions
 * rather than consensus genotypes.
 */

include { BWA_INDEX_SNP           } from '../modules/bwa_index_snp'
include { BWA_MEM_SNP             } from '../modules/bwa_mem_snp'
include { EXTRACT_CODON_FREQS     } from '../modules/extract_codon_freqs'
include { SUMMARIZE_SNP_PROFILING } from '../modules/summarize_snp_profiling'


workflow SNP_PROFILING {

    take:
        ch_reads          // tuple(sample, r1, r2)
        ch_cds_fasta      // path to CDS nucleotide FASTA
        ch_positions_csv  // path to positions of interest CSV

    main:
        // Step 1: Build BWA index (runs once, reused across all samples)
        BWA_INDEX_SNP(ch_cds_fasta)

        // Step 2: Align short reads to CDS reference
        BWA_MEM_SNP(ch_reads, BWA_INDEX_SNP.out.index)

        // Step 3: Extract codon frequencies at target positions
        EXTRACT_CODON_FREQS(BWA_MEM_SNP.out.bam, ch_cds_fasta, ch_positions_csv)

        // Step 4: Generate MultiQC custom content files
        SUMMARIZE_SNP_PROFILING(EXTRACT_CODON_FREQS.out.codon_freqs)

    emit:
        multiqc_snp = SUMMARIZE_SNP_PROFILING.out.generalstats
            .mix(SUMMARIZE_SNP_PROFILING.out.detail)
        codon_freqs = EXTRACT_CODON_FREQS.out.codon_freqs
        versions    = BWA_INDEX_SNP.out.versions
            .mix(BWA_MEM_SNP.out.versions,
                 EXTRACT_CODON_FREQS.out.versions,
                 SUMMARIZE_SNP_PROFILING.out.versions)
}
