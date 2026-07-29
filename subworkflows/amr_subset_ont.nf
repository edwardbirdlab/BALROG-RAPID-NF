/*
 * ONT AMR Subset Detection Subworkflow
 *
 * Read subset approach for long reads:
 *   1. Diamond BLASTX → identify AMR-related reads
 *   2. seqtk subseq → extract matching reads
 *   3. Flye metagenomic assembly → assemble only the AMR subset reads
 *   4. AMRFinderPlus → annotate resistance genes on assembled contigs
 *
 * Handles edge cases: no Diamond hits → empty reads → empty contigs → header-only AMRFinder output
 */

include { DIAMOND_BLASTX_SE } from '../modules/diamond_blastx_se'
include { SEQTK_SUBSEQ_SE  } from '../modules/seqtk_subseq_se'
include { FLYE              } from '../modules/flye'
include { AMRFINDER         } from '../modules/amrfinder'


workflow AMR_SUBSET_ONT {

    take:
        ch_trimmed_reads    // tuple(sample, reads)
        ch_diamond_db       // path to Diamond AMR database (.dmnd)

    main:
        // Step 1: Diamond BLASTX on long reads
        DIAMOND_BLASTX_SE(ch_trimmed_reads, ch_diamond_db)

        // Step 2: Join read IDs back with original reads for extraction
        ch_for_extract = ch_trimmed_reads
            .join(DIAMOND_BLASTX_SE.out.read_ids)

        SEQTK_SUBSEQ_SE(ch_for_extract)

        // Step 3: Flye metagenomic assembly of extracted AMR reads
        FLYE(SEQTK_SUBSEQ_SE.out.subset_reads)

        // Step 4: AMRFinder annotation on assembled contigs
        AMRFINDER(FLYE.out.contigs)

    emit:
        diamond_hits      = DIAMOND_BLASTX_SE.out.hits       // tuple(sample, diamond_tsv)
        subset_reads      = SEQTK_SUBSEQ_SE.out.subset_reads // tuple(sample, reads)
        contigs           = FLYE.out.contigs                  // tuple(sample, contigs_fasta)
        amrfinder_results = AMRFINDER.out.results             // tuple(sample, amrfinder_tsv)
        versions          = DIAMOND_BLASTX_SE.out.versions.mix(SEQTK_SUBSEQ_SE.out.versions, FLYE.out.versions, AMRFINDER.out.versions)
}
