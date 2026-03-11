/*
 * Host Profiling Subworkflow
 *
 * Runs Kraken2 against one or more host genome databases for
 * contamination reporting. Does NOT filter reads.
 *
 * Each sample is classified against each host database (N x M combinations).
 * Outputs: {sample}_{host_name}_k2report.tsv
 */

include { KRAKEN2 as KRAKEN2_HOST } from '../modules/kraken2'


workflow HOST_PROFILING {

    take:
        ch_trimmed_reads    // tuple(sample, r1, r2)
        ch_host_dbs         // tuple(host_name, path_to_db)

    main:
        // Combine every sample with every host database using multiMap
        // to keep the three process inputs synchronized
        ch_combined = ch_trimmed_reads
            .combine(ch_host_dbs)
            .multiMap { sample, r1, r2, host_name, host_db ->
                reads:   tuple(sample, r1, r2)
                db:      host_db
                db_name: host_name
            }

        // Reuse KRAKEN2 module with host db_name → outputs {sample}_{host_name}_k2report.tsv
        KRAKEN2_HOST(
            ch_combined.reads,
            ch_combined.db,
            ch_combined.db_name
        )

    emit:
        host_reports = KRAKEN2_HOST.out.report  // tuple(sample, host_name, k2report)
        versions     = KRAKEN2_HOST.out.versions
}
