/*
 * Read Quality Control Subworkflow
 *
 * Runs FastQC on raw reads, trims with FASTP, optionally trims with BBDuk
 * (for Element Biosciences Aviti reads), then FastQC on final trimmed reads.
 */

include { FASTQC as FASTQC_RAW     } from '../modules/fastqc'
include { FASTP                     } from '../modules/fastp'
include { BBDUK                     } from '../modules/bbduk'
include { FASTQC as FASTQC_TRIMMED } from '../modules/fastqc'


workflow READ_QC {

    take:
        ch_raw_reads       // tuple(sample, r1, r2)
        ch_adapter_fasta   // path to adapter FASTA for BBDuk (value channel)

    main:
        ch_bbduk_log = Channel.empty()

        FASTQC_RAW(ch_raw_reads)
        FASTP(ch_raw_reads)

        if (params.run_bbduk) {
            BBDUK(FASTP.out.trimmed_fastq, ch_adapter_fasta)
            ch_final_reads = BBDUK.out.trimmed
            ch_bbduk_log   = BBDUK.out.log
        } else {
            ch_final_reads = FASTP.out.trimmed_fastq
        }

        FASTQC_TRIMMED(ch_final_reads)

        // Collect versions — include BBDuk when enabled
        ch_versions = FASTQC_RAW.out.versions.mix(FASTP.out.versions, FASTQC_TRIMMED.out.versions)
        if (params.run_bbduk) {
            ch_versions = ch_versions.mix(BBDUK.out.versions)
        }

    emit:
        trimmed_fastq   = ch_final_reads
        fastp_json      = FASTP.out.json
        fastp_html      = FASTP.out.html
        raw_qc_reports  = FASTQC_RAW.out.reports
        trim_qc_reports = FASTQC_TRIMMED.out.reports
        raw_zip         = FASTQC_RAW.out.zip         // zip files for MultiQC
        trim_zip        = FASTQC_TRIMMED.out.zip     // zip files for MultiQC
        bbduk_log       = ch_bbduk_log               // BBDuk log for MultiQC (empty if disabled)
        versions        = ch_versions
}
