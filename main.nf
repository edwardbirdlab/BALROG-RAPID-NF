#!/usr/bin/env nextflow

/*
 * BALROG-RAPID Nextflow Pipeline
 *
 * Metagenomic pathogen detection, host profiling, and AMR detection.
 *
 * Usage:
 *   nextflow run main.nf \
 *     --sample_sheet samplesheet.csv \
 *     --kraken2_db /path/to/kraken2_db \
 *     --sylph_db /path/to/sylph.syldb \
 *     --diamond_db /path/to/amr.dmnd \
 *     --outdir ./results
 *
 * Optional host profiling:
 *   --host_sheet hosts.csv   (CSV with columns: host_name,kraken2_db)
 */

nextflow.enable.dsl = 2

include { BALROG_SHORT_READ } from './workflows/balrog_short_read'


// -------------------------------------------------------------------
// Parameter validation
// -------------------------------------------------------------------

if (!params.sample_sheet) {
    error "ERROR: --sample_sheet is required. Provide a CSV with columns: sample,r1,r2"
}

if (params.run_taxonomy && !params.kraken2_db) {
    error "ERROR: --kraken2_db is required when taxonomy is enabled"
}

if (params.run_taxonomy && !params.sylph_db) {
    error "ERROR: --sylph_db is required when taxonomy is enabled"
}

if (params.run_amr && !params.diamond_db) {
    error "ERROR: --diamond_db is required when AMR detection is enabled"
}

if (params.run_taxonomy && !params.sylph_tax_db) {
    log.warn "WARNING: --sylph_tax_db is not set. Sylph taxonomy profiles (.sylphmpa) will not be generated for MultiQC."
}

if (params.run_bbduk && !params.bbduk_adapters) {
    error "ERROR: --bbduk_adapters is required when BBDuk trimming is enabled (--run_bbduk)"
}

if (params.run_spike_in && !params.spike_in_bt2) {
    error "ERROR: --spike_in_bt2 is required when spike-in removal is enabled (--run_spike_in)"
}


// -------------------------------------------------------------------
// Channel setup
// -------------------------------------------------------------------

// Parse sequencing samplesheet: sample,r1,r2
ch_raw_reads = Channel
    .fromPath(params.sample_sheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        def sample = row.sample
        def r1 = file(row.r1, checkIfExists: true)
        def r2 = file(row.r2, checkIfExists: true)
        tuple(sample, r1, r2)
    }

// Database channels
ch_kraken2_db = params.kraken2_db ? Channel.fromPath(params.kraken2_db, checkIfExists: true).first() : Channel.empty()
ch_sylph_db   = params.sylph_db   ? Channel.fromPath(params.sylph_db,   checkIfExists: true).first() : Channel.empty()
ch_diamond_db   = params.diamond_db   ? Channel.fromPath(params.diamond_db,   checkIfExists: true).first() : Channel.empty()
ch_sylph_tax_db = params.sylph_tax_db ? Channel.fromPath(params.sylph_tax_db, checkIfExists: true).first() : Channel.empty()

// BBDuk adapter channel (optional, for Element Biosciences Aviti runs)
ch_bbduk_adapters = params.run_bbduk && params.bbduk_adapters
    ? Channel.fromPath(params.bbduk_adapters).first()
    : Channel.value([])

// Spike-in Bowtie2 index channel (optional, for T. thermophilus spike-in removal)
ch_spike_in_bt2 = params.run_spike_in && params.spike_in_bt2
    ? Channel.fromPath(params.spike_in_bt2, checkIfExists: true).first()
    : Channel.value([])

// Host database channel (optional)
// CSV format: host_name,kraken2_db
ch_host_dbs = params.host_sheet
    ? Channel
        .fromPath(params.host_sheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.host_name, file(row.kraken2_db, checkIfExists: true))
        }
    : Channel.empty()


// -------------------------------------------------------------------
// Main workflow
// -------------------------------------------------------------------

workflow {

    log.info """
    ======================================
     BALROG-RAPID Pipeline v1.0
    ======================================
     Sample sheet : ${params.sample_sheet}
     Output dir   : ${params.outdir}
     Kraken2 DB   : ${params.kraken2_db ?: 'not set'}
     Sylph DB     : ${params.sylph_db ?: 'not set'}
     Sylph-tax DB : ${params.sylph_tax_db ?: 'not set'}
     Sylph-tax TX : ${params.sylph_tax_taxonomy}
     Diamond DB   : ${params.diamond_db ?: 'not set'}
     Host sheet   : ${params.host_sheet ?: 'not set'}
     Run QC       : ${params.run_qc}
     Run BBDuk    : ${params.run_bbduk}
     BBDuk Adapt  : ${params.run_bbduk ? (params.bbduk_adapters ?: 'not set') : 'N/A'}
     Run Spike-in : ${params.run_spike_in}
     Spike-in Idx : ${params.run_spike_in ? (params.spike_in_bt2 ?: 'not set') : 'N/A'}
     SLURM Acct   : ${params.slurm_account ?: 'not set'}
     Run Taxonomy : ${params.run_taxonomy}
     Run Host     : ${params.run_host_profiling}
     Run AMR      : ${params.run_amr}
     Run MultiQC  : ${params.run_multiqc}
    ======================================
    """.stripIndent()

    BALROG_SHORT_READ(
        ch_raw_reads,
        ch_kraken2_db,
        ch_sylph_db,
        ch_diamond_db,
        ch_host_dbs,
        ch_sylph_tax_db,
        ch_bbduk_adapters,
        ch_spike_in_bt2
    )
}
