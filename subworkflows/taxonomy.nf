/*
 * Taxonomy Profiling Subworkflow
 *
 * Runs Kraken2 (unified database) and Sylph in parallel.
 * Outputs are named for Phase 2 compatibility:
 *   {sample}_k2report.tsv, {sample}_profile.tsv
 */

include { KRAKEN2 } from '../modules/kraken2'
include { SYLPH   } from '../modules/sylph'


workflow TAXONOMY {

    take:
        ch_trimmed_reads    // tuple(sample, r1, r2)
        ch_kraken2_db       // path to unified Kraken2 database directory
        ch_sylph_db         // path to Sylph database file (.syldb)

    main:
        // Kraken2 with empty db_name → outputs {sample}_k2report.tsv (Phase 2 compatible)
        KRAKEN2(ch_trimmed_reads, ch_kraken2_db, "")

        // Sylph profiling → outputs {sample}_profile.tsv
        SYLPH(ch_trimmed_reads, ch_sylph_db)

    emit:
        kraken2_report = KRAKEN2.out.report     // tuple(sample, db_name, k2report)
        kraken2_output = KRAKEN2.out.output     // tuple(sample, db_name, k2out)
        sylph_profile  = SYLPH.out.profile      // tuple(sample, profile)
        versions       = KRAKEN2.out.versions.mix(SYLPH.out.versions)
}
