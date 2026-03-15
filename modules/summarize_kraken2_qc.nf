/*
 * SUMMARIZE_KRAKEN2_QC
 *
 * Parses a Kraken2 report to compute the percentage of reads classified
 * as Bacteria.  Outputs a MultiQC Custom Content TSV that adds a
 * "% Bacterial" column to the General Stats table.
 *
 * Enabled with --custom_qc (requires taxonomy to be running).
 */

process SUMMARIZE_KRAKEN2_QC {

    label 'ultralow'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    input:
        tuple val(sample), val(db_name), path(k2report)

    output:
        path("${sample}_kraken2_qc_mqc.tsv"), emit: generalstats
        path("versions.yml"),                  emit: versions

    script:
    """
    python3 << 'PYTHON_SCRIPT'
import csv

sample = "${sample}"
k2report = "${k2report}"

# Kraken2 report columns (no header):
#   0: percent, 1: clade_reads, 2: taxon_reads, 3: rank, 4: taxid, 5: name
unclassified_reads = 0
root_reads = 0
bacteria_reads = 0

with open(k2report, "r") as fh:
    reader = csv.reader(fh, delimiter="\\t")
    for row in reader:
        if len(row) < 6:
            continue
        rank = row[3].strip()
        taxid = row[4].strip()
        clade = int(row[1].strip())

        if rank == "U" and taxid == "0":
            unclassified_reads = clade
        elif rank == "R" and taxid == "1":
            root_reads = clade
        elif rank == "D" and taxid == "2":
            bacteria_reads = clade

total_reads = unclassified_reads + root_reads
pct_bacterial = (bacteria_reads / total_reads * 100) if total_reads > 0 else 0.0

with open(f"{sample}_kraken2_qc_mqc.tsv", "w") as fh:
    fh.write("# id: 'kraken2_custom_qc'\\n")
    fh.write("# plot_type: 'generalstats'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     - Bacterial_Pct:\\n")
    fh.write("#         title: '% Bacterial'\\n")
    fh.write("#         description: 'Percentage of reads classified as Bacteria by Kraken2'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         max: 100\\n")
    fh.write("#         scale: 'RdYlGn'\\n")
    fh.write("#         suffix: '%'\\n")
    fh.write("#         placement: 1500\\n")
    fh.write("Sample\\tBacterial_Pct\\n")
    fh.write(f"{sample}\\t{pct_bacterial:.2f}\\n")

PYTHON_SCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
