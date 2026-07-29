/*
 * ONT Host Profiling Subworkflow
 *
 * Runs Kraken2 (single-end) against one or more host genome databases
 * for contamination reporting. Does NOT filter reads.
 *
 * Each sample is classified against each host database (N x M combinations).
 * Outputs: {sample}_{host_name}_k2report.tsv
 */

include { KRAKEN2_SE as KRAKEN2_HOST_SE } from '../modules/kraken2_se'


workflow HOST_PROFILING_ONT {

    take:
        ch_trimmed_reads    // tuple(sample, reads)
        ch_host_dbs         // tuple(host_name, path_to_db)

    main:
        // Combine every sample with every host database using multiMap
        // to keep the process inputs synchronized
        ch_combined = ch_trimmed_reads
            .combine(ch_host_dbs)
            .multiMap { sample, reads, host_name, host_db ->
                reads:   tuple(sample, reads)
                db:      host_db
                db_name: host_name
            }

        // Reuse KRAKEN2_SE module with host db_name → outputs {sample}_{host_name}_k2report.tsv
        KRAKEN2_HOST_SE(
            ch_combined.reads,
            ch_combined.db,
            ch_combined.db_name
        )

    emit:
        host_reports = KRAKEN2_HOST_SE.out.report  // tuple(sample, host_name, k2report)
        versions     = KRAKEN2_HOST_SE.out.versions
}
