/*
 * BALROG Short Read Pipeline Workflow
 *
 * Orchestrates all subworkflows:
 *   1. Read QC (FASTP + FastQC)
 *   2. Taxonomy profiling (Kraken2 + Sylph) - parallel
 *   3. Host profiling (Kraken2 x N hosts) - parallel
 *   4. AMR subset detection (Diamond → seqtk → SPAdes → AMRFinder) - parallel
 *   5. Collect software versions from all steps
 *
 * Steps 2-4 all run in parallel after QC completes.
 */

include { READ_QC          } from '../subworkflows/read_qc'
include { TAXONOMY         } from '../subworkflows/taxonomy'
include { HOST_PROFILING   } from '../subworkflows/host_profiling'
include { AMR_SUBSET       } from '../subworkflows/amr_subset'
include { COLLECT_VERSIONS } from '../modules/collect_versions'


workflow BALROG_SHORT_READ {

    take:
        ch_raw_reads    // tuple(sample, r1, r2)
        ch_kraken2_db   // path to unified Kraken2 database
        ch_sylph_db     // path to Sylph database
        ch_diamond_db   // path to Diamond AMR database
        ch_host_dbs     // tuple(host_name, path_to_db) - can be empty channel

    main:
        // Collect all versions.yml files from every process
        ch_versions = Channel.empty()

        // Step 1: Quality control and trimming
        if (params.run_qc) {
            READ_QC(ch_raw_reads)
            ch_reads    = READ_QC.out.trimmed_fastq
            ch_versions = ch_versions.mix(READ_QC.out.versions)
        } else {
            ch_reads = ch_raw_reads
        }

        // Step 2: Taxonomy profiling (runs in parallel with 3 & 4)
        if (params.run_taxonomy) {
            TAXONOMY(ch_reads, ch_kraken2_db, ch_sylph_db)
            ch_versions = ch_versions.mix(TAXONOMY.out.versions)
        }

        // Step 3: Host profiling (runs in parallel with 2 & 4)
        if (params.run_host_profiling) {
            HOST_PROFILING(ch_reads, ch_host_dbs)
            ch_versions = ch_versions.mix(HOST_PROFILING.out.versions)
        }

        // Step 4: AMR detection (runs in parallel with 2 & 3)
        if (params.run_amr) {
            AMR_SUBSET(ch_reads, ch_diamond_db)
            ch_versions = ch_versions.mix(AMR_SUBSET.out.versions)
        }

        // Step 5: Combine all software versions
        COLLECT_VERSIONS(ch_versions.collect())
}
