/*
 * COI Insect ID Subworkflow (short reads only)
 *
 * "The pipeline is: run KMA, filter two columns. Nothing else." -- KMA
 * against the BOLD COI Diptera reference set, thresholded by coverage AND
 * identity (SUMMARIZE_COI). No taxonomy tool, no rank-by-identity scheme:
 * below species, rank is not recoverable from percent identity on COI (see
 * docs/usage.md for the measured evidence).
 *
 * The lineage lookup (processid -> full BOLD lineage) is optional and
 * built once, not per-sample, when --coi_lineage_table is set; otherwise a
 * bundled placeholder is used so SUMMARIZE_COI always has a real file to
 * read regardless.
 */

include { KMA                  } from '../modules/kma'
include { MAKE_LINEAGE_LOOKUP  } from '../modules/make_lineage_lookup'
include { SUMMARIZE_COI        } from '../modules/summarize_coi'


workflow COI_ID {

    take:
        ch_reads          // tuple(sample, r1, r2)
        ch_kma_db         // directory containing the KMA index (value channel)
        ch_lineage_table  // path to the corrected BOLD derep table, or [] (value channel)

    main:
        ch_versions = Channel.empty()

        KMA(ch_reads, ch_kma_db)
        ch_versions = ch_versions.mix(KMA.out.versions)

        if (params.coi_lineage_table) {
            MAKE_LINEAGE_LOOKUP(ch_lineage_table)
            ch_versions = ch_versions.mix(MAKE_LINEAGE_LOOKUP.out.versions)
            ch_lineage_lookup = MAKE_LINEAGE_LOOKUP.out.lookup
        } else {
            ch_lineage_lookup = Channel.fromPath("${projectDir}/assets/no_lineage_lookup.tsv").first()
        }

        SUMMARIZE_COI(KMA.out.res, ch_lineage_lookup)
        ch_versions = ch_versions.mix(SUMMARIZE_COI.out.versions)

        ch_multiqc_coi = SUMMARIZE_COI.out.generalstats.mix(SUMMARIZE_COI.out.detail)

    emit:
        versions    = ch_versions
        multiqc_coi = ch_multiqc_coi
}
