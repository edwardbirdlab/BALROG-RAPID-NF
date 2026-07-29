/*
 * BALROG Long Read (ONT) Pipeline Workflow
 *
 * Orchestrates all ONT subworkflows:
 *   1. Read QC (Porechop_ABI + Filtlong + NanoPlot + FastQC)
 *   2. Taxonomy profiling (Kraken2 SE + Sylph SE) - parallel
 *   3. Host profiling (Kraken2 SE x N hosts) - parallel
 *   4. AMR subset detection (Diamond SE → seqtk SE → Flye → AMRFinder) - parallel
 *   5. Nonpareil coverage estimation (optional)
 *      - When taxonomy enabled: bacterial read extraction → Nonpareil (after step 2)
 *      - When taxonomy disabled: runs on full reads (parallel with 2-4)
 *
 * Steps 2-4 run in parallel. Step 5 runs after step 2 when taxonomy is enabled.
 * No spike-in removal step (short-read specific).
 */

include { READ_QC_ONT        } from '../subworkflows/read_qc_ont'
include { TAXONOMY_ONT       } from '../subworkflows/taxonomy_ont'
include { HOST_PROFILING_ONT } from '../subworkflows/host_profiling_ont'
include { AMR_SUBSET_ONT     } from '../subworkflows/amr_subset_ont'
include { EXTRACT_BACTERIAL_READS_SE } from '../modules/extract_bacterial_reads_se'
include { NONPAREIL_SE               } from '../modules/nonpareil_se'
include { SUMMARIZE_AMRFINDER  as SUMMARIZE_AMRFINDER_ONT  } from '../modules/summarize_amrfinder'
include { SUMMARIZE_KRAKEN2_QC as SUMMARIZE_KRAKEN2_QC_ONT } from '../modules/summarize_kraken2_qc'
include { SNP_PROFILING_ONT } from '../subworkflows/snp_profiling_ont'


workflow BALROG_LONG_READ {

    take:
        ch_raw_reads      // tuple(sample, reads)
        ch_kraken2_db     // path to unified Kraken2 database
        ch_sylph_db       // path to Sylph database
        ch_diamond_db     // path to Diamond AMR database
        ch_host_dbs        // tuple(host_name, path_to_db) - can be empty channel
        ch_sylph_tax_db    // path to pre-downloaded sylph-tax taxonomy DB directory
        ch_snp_cds_fasta   // path to CDS FASTA for SNP profiling (value channel)
        ch_snp_positions_csv // path to positions CSV for SNP profiling (value channel)

    main:
        // Collect all versions.yml files from every process
        ch_versions = Channel.empty()

        // Initialize MultiQC collection channels (empty defaults for disabled steps)
        ch_multiqc_fastqc_raw  = Channel.empty()
        ch_multiqc_fastqc_trim = Channel.empty()
        ch_multiqc_nanoplot    = Channel.empty()
        ch_multiqc_porechop    = Channel.empty()
        ch_multiqc_filtlong    = Channel.empty()
        ch_multiqc_k2_taxonomy = Channel.empty()
        ch_multiqc_k2_host     = Channel.empty()
        ch_multiqc_sylph       = Channel.empty()
        ch_multiqc_amrfinder   = Channel.empty()
        ch_multiqc_nonpareil   = Channel.empty()
        ch_multiqc_custom_qc   = Channel.empty()
        ch_multiqc_snp         = Channel.empty()

        // Step 1: Quality control and trimming (ONT)
        if (params.run_qc) {
            READ_QC_ONT(ch_raw_reads)
            ch_reads    = READ_QC_ONT.out.trimmed_reads
            ch_versions = ch_versions.mix(READ_QC_ONT.out.versions)

            // Collect QC outputs for MultiQC
            ch_multiqc_fastqc_raw  = READ_QC_ONT.out.raw_zip
            ch_multiqc_fastqc_trim = READ_QC_ONT.out.trim_zip
            ch_multiqc_nanoplot    = READ_QC_ONT.out.nanoplot_stats
            ch_multiqc_porechop    = READ_QC_ONT.out.porechop_log
            ch_multiqc_filtlong    = READ_QC_ONT.out.filtlong_log
        } else {
            ch_reads = ch_raw_reads
        }

        // Step 2: Taxonomy profiling (runs in parallel with 3 & 4)
        if (params.run_taxonomy) {
            TAXONOMY_ONT(ch_reads, ch_kraken2_db, ch_sylph_db, ch_sylph_tax_db)
            ch_versions = ch_versions.mix(TAXONOMY_ONT.out.versions)

            // Extract file paths from tuples for MultiQC
            ch_multiqc_k2_taxonomy = TAXONOMY_ONT.out.kraken2_report.map { it[-1] }
            ch_multiqc_sylph       = TAXONOMY_ONT.out.sylph_tax_mpa

            // Custom QC: bacterial read percentage for MultiQC General Stats
            if (params.custom_qc) {
                SUMMARIZE_KRAKEN2_QC_ONT(TAXONOMY_ONT.out.kraken2_report)
                ch_versions = ch_versions.mix(SUMMARIZE_KRAKEN2_QC_ONT.out.versions)
                ch_multiqc_custom_qc = SUMMARIZE_KRAKEN2_QC_ONT.out.generalstats
            }
        }

        // Step 3: Host profiling (runs in parallel with 2 & 4)
        if (params.run_host_profiling) {
            HOST_PROFILING_ONT(ch_reads, ch_host_dbs)
            ch_versions = ch_versions.mix(HOST_PROFILING_ONT.out.versions)

            // Extract file paths from tuples for MultiQC
            ch_multiqc_k2_host = HOST_PROFILING_ONT.out.host_reports.map { it[-1] }
        }

        // Step 4: AMR detection (runs in parallel with 2 & 3)
        if (params.run_amr) {
            AMR_SUBSET_ONT(ch_reads, ch_diamond_db)
            ch_versions = ch_versions.mix(AMR_SUBSET_ONT.out.versions)

            // Summarize AMRFinder results for MultiQC custom content
            SUMMARIZE_AMRFINDER_ONT(AMR_SUBSET_ONT.out.amrfinder_results)
            ch_versions = ch_versions.mix(SUMMARIZE_AMRFINDER_ONT.out.versions)
            ch_multiqc_amrfinder = SUMMARIZE_AMRFINDER_ONT.out.generalstats
                .mix(SUMMARIZE_AMRFINDER_ONT.out.classes)
                .mix(SUMMARIZE_AMRFINDER_ONT.out.detail)
        }

        // Step 5: Nonpareil coverage estimation
        //   When taxonomy enabled: extract bacterial reads first (removes host bias)
        //   When taxonomy disabled: run on full reads (original behavior)
        if (params.run_nonpareil) {
            if (params.run_taxonomy) {
                // Strip db_name from Kraken2 outputs, join with original reads
                ch_k2_output_stripped = TAXONOMY_ONT.out.kraken2_output
                    .map { sample, db_name, k2out -> tuple(sample, k2out) }
                ch_k2_report_stripped = TAXONOMY_ONT.out.kraken2_report
                    .map { sample, db_name, k2report -> tuple(sample, k2report) }

                ch_for_extraction = ch_reads
                    .join(ch_k2_output_stripped)
                    .join(ch_k2_report_stripped)
                    // Result: tuple(sample, reads, k2out, k2report)

                EXTRACT_BACTERIAL_READS_SE(ch_for_extraction)
                ch_versions = ch_versions.mix(EXTRACT_BACTERIAL_READS_SE.out.versions)

                // Feed bacterial reads to Nonpareil
                NONPAREIL_SE(EXTRACT_BACTERIAL_READS_SE.out.bacterial_reads)
            } else {
                // No taxonomy data available — run Nonpareil on full reads
                NONPAREIL_SE(ch_reads)
            }
            ch_versions = ch_versions.mix(NONPAREIL_SE.out.versions)
            ch_multiqc_nonpareil = NONPAREIL_SE.out.json
        }

        // Step 6: Targeted SNP/AA variant profiling (runs in parallel with 2-4)
        if (params.run_snp_profiling) {
            SNP_PROFILING_ONT(ch_reads, ch_snp_cds_fasta, ch_snp_positions_csv)
            ch_versions = ch_versions.mix(SNP_PROFILING_ONT.out.versions)
            ch_multiqc_snp = SNP_PROFILING_ONT.out.multiqc_snp
        }

    emit:
        versions             = ch_versions
        multiqc_fastqc_raw   = ch_multiqc_fastqc_raw
        multiqc_fastqc_trim  = ch_multiqc_fastqc_trim
        multiqc_nanoplot     = ch_multiqc_nanoplot
        multiqc_porechop     = ch_multiqc_porechop
        multiqc_filtlong     = ch_multiqc_filtlong
        multiqc_k2_taxonomy  = ch_multiqc_k2_taxonomy
        multiqc_k2_host      = ch_multiqc_k2_host
        multiqc_sylph        = ch_multiqc_sylph
        multiqc_amrfinder    = ch_multiqc_amrfinder
        multiqc_nonpareil    = ch_multiqc_nonpareil
        multiqc_custom_qc    = ch_multiqc_custom_qc
        multiqc_snp          = ch_multiqc_snp
}
