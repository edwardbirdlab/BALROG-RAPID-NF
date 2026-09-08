#!/usr/bin/env nextflow

/*
 * BALROG-RAPID Nextflow Pipeline
 *
 * Unified multi-library metagenomic pipeline supporting both
 * Illumina/Element short reads and Oxford Nanopore long reads.
 *
 * Usage:
 *   nextflow run main.nf \
 *     -profile docker \
 *     --sample_sheet samplesheet.csv \
 *     --kraken2_db /path/to/kraken2_db \
 *     --sylph_db /path/to/sylph.syldb \
 *     --diamond_db /path/to/amr.dmnd \
 *     --outdir ./results
 *
 * Samplesheet format (CSV):
 *   sample,library,molecule,r1,r2,reads
 *   S1,Short,gDNA,/path/S1_R1.fq.gz,/path/S1_R2.fq.gz,
 *   S2,Long,gDNA,,,/path/S2.fq.gz
 *
 * Optional host profiling:
 *   --host_sheet hosts.csv   (CSV with columns: host_name,kraken2_db)
 */

nextflow.enable.dsl = 2

include { BALROG_SHORT_READ } from './workflows/balrog_short_read'
include { BALROG_LONG_READ  } from './workflows/balrog_long_read'
include { COLLECT_VERSIONS  } from './modules/collect_versions'
include { MULTIQC           } from './modules/multiqc'
include { CAT_FASTQ_PE      } from './modules/cat_fastq_pe'
include { CAT_FASTQ_SE      } from './modules/cat_fastq_se'


workflow {

    // -------------------------------------------------------------------
    // Parameter validation
    // -------------------------------------------------------------------

    if (!params.sample_sheet) {
        error "ERROR: --sample_sheet is required. Provide a CSV with columns: sample,library,molecule,r1,r2,reads"
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

    if (params.custom_qc && !params.run_taxonomy) {
        log.warn "WARNING: --custom_qc has no effect without --run_taxonomy enabled"
    }

    if (params.run_snp_profiling && !params.snp_cds_fasta) {
        error "ERROR: --snp_cds_fasta is required when --run_snp_profiling is enabled"
    }
    if (params.run_snp_profiling && !params.snp_positions_csv) {
        error "ERROR: --snp_positions_csv is required when --run_snp_profiling is enabled"
    }

    if (params.run_coi_id && !params.kma_db) {
        error "ERROR: --kma_db is required when --run_coi_id is enabled"
    }


    // -------------------------------------------------------------------
    // Channel setup
    // -------------------------------------------------------------------

    // Parse unified samplesheet: sample,library,molecule,r1,r2,reads
    // Supports re-sequencing: multiple rows per sample are concatenated.
    ch_all_samples = Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true)

    // Split by library type, group by sample for re-sequencing concatenation,
    // and validate that all rows sharing a sample name have consistent library+molecule.
    ch_short_grouped = ch_all_samples
        .filter { it.library == 'Short' }
        .map { row ->
            tuple(row.sample, row.library, row.molecule,
                  file(row.r1, checkIfExists: true), file(row.r2, checkIfExists: true))
        }
        .groupTuple(by: 0)
        .map { sample, libs, mols, r1s, r2s ->
            if (libs.unique().size() > 1 || mols.unique().size() > 1) {
                error "ERROR: Sample '${sample}' has inconsistent library/molecule across " +
                      "re-sequencing rows: library=${libs.unique()}, molecule=${mols.unique()}. " +
                      "All rows for the same sample must share identical library and molecule values."
            }
            tuple(sample, r1s, r2s)
        }

    ch_long_grouped = ch_all_samples
        .filter { it.library == 'Long' }
        .map { row ->
            tuple(row.sample, row.library, row.molecule,
                  file(row.reads, checkIfExists: true))
        }
        .groupTuple(by: 0)
        .map { sample, libs, mols, reads ->
            if (libs.unique().size() > 1 || mols.unique().size() > 1) {
                error "ERROR: Sample '${sample}' has inconsistent library/molecule across " +
                      "re-sequencing rows: library=${libs.unique()}, molecule=${mols.unique()}. " +
                      "All rows for the same sample must share identical library and molecule values."
            }
            tuple(sample, reads)
        }

    // Branch: single-file groups pass through, multi-file groups go to CAT_FASTQ
    ch_short_grouped
        .branch {
            single: it[1].size() == 1
            multi:  it[1].size() > 1
        }
        .set { ch_short_branched }

    ch_long_grouped
        .branch {
            single: it[1].size() == 1
            multi:  it[1].size() > 1
        }
        .set { ch_long_branched }

    // Single-file groups: unwrap from list
    ch_short_single = ch_short_branched.single
        .map { sample, r1s, r2s -> tuple(sample, r1s[0], r2s[0]) }
    ch_long_single = ch_long_branched.single
        .map { sample, reads -> tuple(sample, reads[0]) }

    // Database channels
    ch_kraken2_db   = params.kraken2_db   ? Channel.fromPath(params.kraken2_db,   checkIfExists: true).first() : Channel.empty()
    ch_sylph_db     = params.sylph_db     ? Channel.fromPath(params.sylph_db,     checkIfExists: true).first() : Channel.empty()
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

    // SNP profiling channels (optional)
    ch_snp_cds_fasta = params.run_snp_profiling && params.snp_cds_fasta
        ? Channel.fromPath(params.snp_cds_fasta, checkIfExists: true).first()
        : Channel.value([])
    ch_snp_positions_csv = params.run_snp_profiling && params.snp_positions_csv
        ? Channel.fromPath(params.snp_positions_csv, checkIfExists: true).first()
        : Channel.value([])

    // COI insect ID channels (optional, short-read only)
    ch_kma_db = params.run_coi_id && params.kma_db
        ? Channel.fromPath(params.kma_db, checkIfExists: true).first()
        : Channel.value([])
    // MAKE_LINEAGE_LOOKUP reads this directly only when --coi_lineage_table is
    // set (see subworkflows/coi_id.nf); otherwise it's never consumed.
    ch_coi_lineage_table = params.coi_lineage_table
        ? Channel.fromPath(params.coi_lineage_table, checkIfExists: true).first()
        : Channel.value([])


    // -------------------------------------------------------------------
    // Main workflow
    // -------------------------------------------------------------------

    // Multi-file groups: concatenate
    CAT_FASTQ_PE(ch_short_branched.multi)
    CAT_FASTQ_SE(ch_long_branched.multi)

    // Merge into final channels
    ch_short_reads = ch_short_single.mix(CAT_FASTQ_PE.out.reads)
    ch_long_reads  = ch_long_single.mix(CAT_FASTQ_SE.out.reads)

    log.info """
    ======================================
     BALROG-RAPID Pipeline v2.0
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
     Run Nonpar.  : ${params.run_nonpareil}
     Run SNP Prof : ${params.run_snp_profiling}
     SNP CDS FA   : ${params.run_snp_profiling ? (params.snp_cds_fasta ?: 'not set') : 'N/A'}
     SNP Pos CSV  : ${params.run_snp_profiling ? (params.snp_positions_csv ?: 'not set') : 'N/A'}
     SNP Min BQ   : ${params.run_snp_profiling ? params.snp_min_base_quality : 'N/A'}
     SNP Min MAPQ : ${params.run_snp_profiling ? params.snp_min_mapq : 'N/A'}
     Custom QC    : ${params.custom_qc}
     Run MultiQC  : ${params.run_multiqc}
     Filtlong len : ${params.filtlong_min_length}
     Filtlong pct : ${params.filtlong_keep_percent}
     Flye mode    : ${params.flye_read_type}
     Run COI ID   : ${params.run_coi_id}
     KMA DB       : ${params.run_coi_id ? (params.kma_db ?: 'not set') : 'N/A'}
     COI Lineage  : ${params.coi_lineage_table ?: 'not set (BIN/taxon only, no lineage join)'}
    ======================================
    """.stripIndent()

    // -------------------------------------------------------------------
    // Run short-read workflow (if any Short library samples exist)
    // -------------------------------------------------------------------
    BALROG_SHORT_READ(
        ch_short_reads,
        ch_kraken2_db,
        ch_sylph_db,
        ch_diamond_db,
        ch_host_dbs,
        ch_sylph_tax_db,
        ch_bbduk_adapters,
        ch_spike_in_bt2,
        ch_snp_cds_fasta,
        ch_snp_positions_csv,
        ch_kma_db,
        ch_coi_lineage_table
    )

    // -------------------------------------------------------------------
    // Run long-read workflow (if any Long library samples exist)
    // -------------------------------------------------------------------
    BALROG_LONG_READ(
        ch_long_reads,
        ch_kraken2_db,
        ch_sylph_db,
        ch_diamond_db,
        ch_host_dbs,
        ch_sylph_tax_db,
        ch_snp_cds_fasta,
        ch_snp_positions_csv
    )

    // -------------------------------------------------------------------
    // Merge versions from both workflows and collect
    // -------------------------------------------------------------------
    ch_all_versions = BALROG_SHORT_READ.out.versions
        .mix(BALROG_LONG_READ.out.versions)
        .mix(CAT_FASTQ_PE.out.versions.ifEmpty([]))
        .mix(CAT_FASTQ_SE.out.versions.ifEmpty([]))

    COLLECT_VERSIONS(ch_all_versions.collect())

    // -------------------------------------------------------------------
    // MultiQC: merge all outputs from both workflows into unified report
    // -------------------------------------------------------------------
    if (params.run_multiqc) {
        ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)

        MULTIQC(
            // FastQC (merged short + long)
            BALROG_SHORT_READ.out.multiqc_fastqc_raw
                .mix(BALROG_LONG_READ.out.multiqc_fastqc_raw)
                .collect().ifEmpty([]),
            // FASTP (short-read only)
            BALROG_SHORT_READ.out.multiqc_fastp
                .collect().ifEmpty([]),
            // BBDuk (short-read only)
            BALROG_SHORT_READ.out.multiqc_bbduk
                .collect().ifEmpty([]),
            // FastQC trimmed (merged short + long)
            BALROG_SHORT_READ.out.multiqc_fastqc_trim
                .mix(BALROG_LONG_READ.out.multiqc_fastqc_trim)
                .collect().ifEmpty([]),
            // Spike-in (short-read only)
            BALROG_SHORT_READ.out.multiqc_spike_in
                .collect().ifEmpty([]),
            // Kraken2 taxonomy (merged)
            BALROG_SHORT_READ.out.multiqc_k2_taxonomy
                .mix(BALROG_LONG_READ.out.multiqc_k2_taxonomy)
                .collect().ifEmpty([]),
            // Kraken2 host (merged)
            BALROG_SHORT_READ.out.multiqc_k2_host
                .mix(BALROG_LONG_READ.out.multiqc_k2_host)
                .collect().ifEmpty([]),
            // Sylph (merged)
            BALROG_SHORT_READ.out.multiqc_sylph
                .mix(BALROG_LONG_READ.out.multiqc_sylph)
                .collect().ifEmpty([]),
            // AMRFinder (merged)
            BALROG_SHORT_READ.out.multiqc_amrfinder
                .mix(BALROG_LONG_READ.out.multiqc_amrfinder)
                .collect().ifEmpty([]),
            // Nonpareil (merged)
            BALROG_SHORT_READ.out.multiqc_nonpareil
                .mix(BALROG_LONG_READ.out.multiqc_nonpareil)
                .collect().ifEmpty([]),
            // Custom QC (merged)
            BALROG_SHORT_READ.out.multiqc_custom_qc
                .mix(BALROG_LONG_READ.out.multiqc_custom_qc)
                .collect().ifEmpty([]),
            // NanoPlot (long-read only)
            BALROG_LONG_READ.out.multiqc_nanoplot
                .collect().ifEmpty([]),
            // Porechop (long-read only)
            BALROG_LONG_READ.out.multiqc_porechop
                .collect().ifEmpty([]),
            // Filtlong (long-read only)
            BALROG_LONG_READ.out.multiqc_filtlong
                .collect().ifEmpty([]),
            // SNP Profiling (merged short + long)
            BALROG_SHORT_READ.out.multiqc_snp
                .mix(BALROG_LONG_READ.out.multiqc_snp)
                .collect().ifEmpty([]),
            // COI Insect ID (short-read only)
            BALROG_SHORT_READ.out.multiqc_coi
                .collect().ifEmpty([]),
            // Config and versions
            ch_multiqc_config.first(),
            COLLECT_VERSIONS.out.combined_versions
        )
    }
}
