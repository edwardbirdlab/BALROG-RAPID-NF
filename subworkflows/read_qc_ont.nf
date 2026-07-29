/*
 * ONT Read Quality Control Subworkflow
 *
 * Runs FastQC on raw reads, trims adapters with Porechop_ABI,
 * filters by quality/length with Filtlong, then runs NanoPlot
 * and FastQC on the final trimmed reads.
 */

include { FASTQC_SE as FASTQC_RAW_SE     } from '../modules/fastqc_se'
include { PORECHOP_ABI                    } from '../modules/porechop_abi'
include { FILTLONG                        } from '../modules/filtlong'
include { NANOPLOT                        } from '../modules/nanoplot'
include { FASTQC_SE as FASTQC_TRIMMED_SE } from '../modules/fastqc_se'


workflow READ_QC_ONT {

    take:
        ch_raw_reads    // tuple(sample, reads)

    main:
        // FastQC on raw reads
        FASTQC_RAW_SE(ch_raw_reads)

        // Adapter trimming with Porechop_ABI (ab initio discovery)
        PORECHOP_ABI(ch_raw_reads)

        // Quality and length filtering
        FILTLONG(PORECHOP_ABI.out.trimmed)

        // QC on final trimmed reads
        NANOPLOT(FILTLONG.out.filtered)
        FASTQC_TRIMMED_SE(FILTLONG.out.filtered)

        // Collect versions
        ch_versions = FASTQC_RAW_SE.out.versions
            .mix(PORECHOP_ABI.out.versions,
                 FILTLONG.out.versions,
                 NANOPLOT.out.versions,
                 FASTQC_TRIMMED_SE.out.versions)

    emit:
        trimmed_reads   = FILTLONG.out.filtered
        raw_zip         = FASTQC_RAW_SE.out.zip
        trim_zip        = FASTQC_TRIMMED_SE.out.zip
        nanoplot_stats  = NANOPLOT.out.stats
        porechop_log    = PORECHOP_ABI.out.log
        filtlong_log    = FILTLONG.out.log
        versions        = ch_versions
}
