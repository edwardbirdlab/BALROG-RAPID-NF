/*
 * ONT Taxonomy Profiling Subworkflow
 *
 * Runs Kraken2 (single-end) and Sylph (single-end) in parallel.
 * Optionally runs sylph-tax to convert Sylph output into
 * taxonomic profiles (.sylphmpa) for MultiQC.
 *
 * Outputs are named for Phase 2 compatibility:
 *   {sample}_k2report.tsv, {sample}_profile.tsv
 */

include { KRAKEN2_SE } from '../modules/kraken2_se'
include { SYLPH_SE   } from '../modules/sylph_se'
include { SYLPH_TAX  } from '../modules/sylph_tax'


workflow TAXONOMY_ONT {

    take:
        ch_trimmed_reads    // tuple(sample, reads)
        ch_kraken2_db       // path to unified Kraken2 database directory
        ch_sylph_db         // path to Sylph database file (.syldb)
        ch_sylph_tax_db     // path to pre-downloaded sylph-tax taxonomy DB directory

    main:
        // Kraken2 with empty db_name → outputs {sample}_k2report.tsv (Phase 2 compatible)
        KRAKEN2_SE(ch_trimmed_reads, ch_kraken2_db, "")

        // Sylph profiling → outputs {sample}_profile.tsv
        SYLPH_SE(ch_trimmed_reads, ch_sylph_db)

        // Convert sylph profile → .sylphmpa for MultiQC (uses sylph-tax)
        SYLPH_TAX(SYLPH_SE.out.profile, ch_sylph_tax_db, params.sylph_tax_taxonomy)

    emit:
        kraken2_report = KRAKEN2_SE.out.report      // tuple(sample, db_name, k2report)
        kraken2_output = KRAKEN2_SE.out.output       // tuple(sample, db_name, k2out)
        sylph_profile  = SYLPH_SE.out.profile        // tuple(sample, profile)
        sylph_tax_mpa  = SYLPH_TAX.out.sylphmpa      // path to .sylphmpa files
        versions       = KRAKEN2_SE.out.versions.mix(SYLPH_SE.out.versions, SYLPH_TAX.out.versions)
}
