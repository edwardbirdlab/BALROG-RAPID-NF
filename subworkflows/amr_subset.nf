/*
 * AMR Subset Detection Subworkflow
 *
 * Read subset approach:
 *   1. Diamond BLASTX on R1 reads → identify AMR-related reads
 *   2. seqtk subseq → extract matching read pairs from R1 and R2
 *   3. SPAdes micro-assembly → assemble only the AMR subset reads
 *   4. AMRFinderPlus → annotate resistance genes on assembled contigs
 *
 * Handles edge cases: no Diamond hits → empty reads → empty contigs → header-only AMRFinder output
 */

include { DIAMOND_BLASTX } from '../modules/diamond_blastx'
include { SEQTK_SUBSEQ  } from '../modules/seqtk_subseq'
include { SPADES         } from '../modules/spades'
include { AMRFINDER      } from '../modules/amrfinder'


workflow AMR_SUBSET {

    take:
        ch_trimmed_reads    // tuple(sample, r1, r2)
        ch_diamond_db       // path to Diamond AMR database (.dmnd)

    main:
        // Step 1: Diamond BLASTX on R1 reads
        DIAMOND_BLASTX(ch_trimmed_reads, ch_diamond_db)

        // Step 2: Join read IDs back with original reads for extraction
        ch_for_extract = ch_trimmed_reads
            .join(DIAMOND_BLASTX.out.read_ids)

        SEQTK_SUBSEQ(ch_for_extract)

        // Step 3: Micro-assembly of extracted AMR reads
        SPADES(SEQTK_SUBSEQ.out.subset_reads)

        // Step 4: AMRFinder annotation on assembled contigs
        AMRFINDER(SPADES.out.contigs)

    emit:
        diamond_hits      = DIAMOND_BLASTX.out.hits        // tuple(sample, diamond_tsv)
        subset_reads      = SEQTK_SUBSEQ.out.subset_reads // tuple(sample, r1, r2)
        contigs           = SPADES.out.contigs             // tuple(sample, contigs_fasta)
        amrfinder_results = AMRFINDER.out.results          // tuple(sample, amrfinder_tsv)
        versions          = DIAMOND_BLASTX.out.versions.mix(SEQTK_SUBSEQ.out.versions, SPADES.out.versions, AMRFINDER.out.versions)
}
