/*
 * Read Quality Control Subworkflow
 *
 * Runs FastQC on raw reads, trims with FASTP, then FastQC on trimmed reads.
 */

include { FASTQC as FASTQC_RAW     } from '../modules/fastqc'
include { FASTP                     } from '../modules/fastp'
include { FASTQC as FASTQC_TRIMMED } from '../modules/fastqc'


workflow READ_QC {

    take:
        ch_raw_reads    // tuple(sample, r1, r2)

    main:
        FASTQC_RAW(ch_raw_reads)
        FASTP(ch_raw_reads)
        FASTQC_TRIMMED(FASTP.out.trimmed_fastq)

    emit:
        trimmed_fastq   = FASTP.out.trimmed_fastq    // tuple(sample, r1_trimmed, r2_trimmed)
        fastp_json      = FASTP.out.json
        fastp_html      = FASTP.out.html
        raw_qc_reports  = FASTQC_RAW.out.reports
        trim_qc_reports = FASTQC_TRIMMED.out.reports
        versions        = FASTQC_RAW.out.versions.mix(FASTP.out.versions, FASTQC_TRIMMED.out.versions)
}
