/*
 * BALROG Short Read Pipeline Workflow
 *
 * Orchestrates all subworkflows:
 *   1. Read QC (FASTP + optional BBDuk + FastQC)
 *   1b. Optional spike-in removal (Bowtie2 vs T. thermophilus)
 *   2. Taxonomy profiling (Kraken2 + Sylph) - parallel
 *   3. Host profiling (Kraken2 x N hosts) - parallel
 *   4. AMR subset detection (Diamond → seqtk → SPAdes → AMRFinder) - parallel
 *   5. Collect software versions from all steps
 *   6. MultiQC aggregation report
 *
 * Steps 2-4 all run in parallel after QC (and optional spike-in removal) completes.
 */

include { READ_QC          } from '../subworkflows/read_qc'
include { TAXONOMY         } from '../subworkflows/taxonomy'
include { HOST_PROFILING   } from '../subworkflows/host_profiling'
include { AMR_SUBSET       } from '../subworkflows/amr_subset'
include { SPIKE_IN_REMOVAL      } from '../modules/spike_in_removal'
include { COLLECT_VERSIONS      } from '../modules/collect_versions'
include { SUMMARIZE_AMRFINDER   } from '../modules/summarize_amrfinder'
include { MULTIQC               } from '../modules/multiqc'


workflow BALROG_SHORT_READ {

    take:
        ch_raw_reads      // tuple(sample, r1, r2)
        ch_kraken2_db     // path to unified Kraken2 database
        ch_sylph_db       // path to Sylph database
        ch_diamond_db     // path to Diamond AMR database
        ch_host_dbs       // tuple(host_name, path_to_db) - can be empty channel
        ch_sylph_tax_db   // path to pre-downloaded sylph-tax taxonomy DB directory
        ch_bbduk_adapters // path to adapter FASTA for BBDuk (value channel)
        ch_spike_in_bt2   // path to pre-built Bowtie2 index directory (value channel)

    main:
        // Collect all versions.yml files from every process
        ch_versions = Channel.empty()

        // Initialize MultiQC collection channels (empty defaults for disabled steps)
        ch_multiqc_fastqc_raw  = Channel.empty()
        ch_multiqc_fastp       = Channel.empty()
        ch_multiqc_bbduk       = Channel.empty()
        ch_multiqc_fastqc_trim = Channel.empty()
        ch_multiqc_spike_in    = Channel.empty()
        ch_multiqc_k2_taxonomy = Channel.empty()
        ch_multiqc_k2_host     = Channel.empty()
        ch_multiqc_sylph       = Channel.empty()
        ch_multiqc_amrfinder   = Channel.empty()

        // Step 1: Quality control and trimming
        if (params.run_qc) {
            READ_QC(ch_raw_reads, ch_bbduk_adapters)
            ch_reads    = READ_QC.out.trimmed_fastq
            ch_versions = ch_versions.mix(READ_QC.out.versions)

            // Collect QC outputs for MultiQC
            ch_multiqc_fastqc_raw  = READ_QC.out.raw_zip
            ch_multiqc_fastp       = READ_QC.out.fastp_json
            ch_multiqc_bbduk       = READ_QC.out.bbduk_log
            ch_multiqc_fastqc_trim = READ_QC.out.trim_zip
        } else {
            ch_reads = ch_raw_reads
        }

        // Step 1b: Spike-in removal (optional, after QC, before downstream analysis)
        if (params.run_spike_in) {
            SPIKE_IN_REMOVAL(ch_reads, ch_spike_in_bt2)
            ch_reads    = SPIKE_IN_REMOVAL.out.cleaned_reads
            ch_versions = ch_versions.mix(SPIKE_IN_REMOVAL.out.versions)
            ch_multiqc_spike_in = SPIKE_IN_REMOVAL.out.log
        }

        // Step 2: Taxonomy profiling (runs in parallel with 3 & 4)
        if (params.run_taxonomy) {
            TAXONOMY(ch_reads, ch_kraken2_db, ch_sylph_db, ch_sylph_tax_db)
            ch_versions = ch_versions.mix(TAXONOMY.out.versions)

            // Extract file paths from tuples for MultiQC
            ch_multiqc_k2_taxonomy = TAXONOMY.out.kraken2_report.map { it[-1] }
            ch_multiqc_sylph       = TAXONOMY.out.sylph_tax_mpa
        }

        // Step 3: Host profiling (runs in parallel with 2 & 4)
        if (params.run_host_profiling) {
            HOST_PROFILING(ch_reads, ch_host_dbs)
            ch_versions = ch_versions.mix(HOST_PROFILING.out.versions)

            // Extract file paths from tuples for MultiQC
            ch_multiqc_k2_host = HOST_PROFILING.out.host_reports.map { it[-1] }
        }

        // Step 4: AMR detection (runs in parallel with 2 & 3)
        if (params.run_amr) {
            AMR_SUBSET(ch_reads, ch_diamond_db)
            ch_versions = ch_versions.mix(AMR_SUBSET.out.versions)

            // Summarize AMRFinder results for MultiQC custom content
            SUMMARIZE_AMRFINDER(AMR_SUBSET.out.amrfinder_results)
            ch_versions = ch_versions.mix(SUMMARIZE_AMRFINDER.out.versions)
            ch_multiqc_amrfinder = SUMMARIZE_AMRFINDER.out.generalstats
                .mix(SUMMARIZE_AMRFINDER.out.classes)
                .mix(SUMMARIZE_AMRFINDER.out.detail)
        }

        // Step 5: Combine all software versions
        COLLECT_VERSIONS(ch_versions.collect())

        // Step 6: MultiQC aggregation report
        if (params.run_multiqc) {
            ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)

            MULTIQC(
                ch_multiqc_fastqc_raw.collect().ifEmpty([]),
                ch_multiqc_fastp.collect().ifEmpty([]),
                ch_multiqc_bbduk.collect().ifEmpty([]),
                ch_multiqc_fastqc_trim.collect().ifEmpty([]),
                ch_multiqc_spike_in.collect().ifEmpty([]),
                ch_multiqc_k2_taxonomy.collect().ifEmpty([]),
                ch_multiqc_k2_host.collect().ifEmpty([]),
                ch_multiqc_sylph.collect().ifEmpty([]),
                ch_multiqc_amrfinder.collect().ifEmpty([]),
                ch_multiqc_config.first(),
                COLLECT_VERSIONS.out.combined_versions
            )
        }
}
