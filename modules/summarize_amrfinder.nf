/*
 * SUMMARIZE_AMRFINDER
 *
 * Reads an AMRFinderPlus TSV and produces three MultiQC Custom Content files:
 *   1. General stats  (_amrfinder_generalstats_mqc.tsv)
 *   2. AMR class bar   (_amrfinder_classes_mqc.tsv)
 *   3. Detail table    (_amrfinder_detail_mqc.tsv)
 *
 * These _mqc.tsv files are auto-detected by MultiQC's Custom Content module.
 */

process SUMMARIZE_AMRFINDER {

    label 'ultralow'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    input:
        tuple val(sample), path(amrfinder_tsv)

    output:
        path("${sample}_amrfinder_generalstats_mqc.tsv"), emit: generalstats
        path("${sample}_amrfinder_classes_mqc.tsv"),      emit: classes
        path("${sample}_amrfinder_detail_mqc.tsv"),       emit: detail
        path("versions.yml"),                              emit: versions

    script:
    """
    python3 << 'PYTHON_SCRIPT'
import csv
import sys
from collections import Counter

sample = "${sample}"
tsv_path = "${amrfinder_tsv}"

# --- Read AMRFinder TSV ---
rows = []
with open(tsv_path, "r") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        rows.append(row)

# --- Categorize genes by Type ---
type_counts = Counter()
amr_classes = Counter()
methods = Counter()
coverages = []
identities = []

for row in rows:
    gene_type = row.get("Type", "Unknown")
    type_counts[gene_type] += 1

    # Count AMR drug classes
    if gene_type == "AMR":
        drug_class = row.get("Class", "Unknown")
        if drug_class:
            amr_classes[drug_class] += 1

    # Method breakdown
    method = row.get("Method", "Unknown")
    if method:
        methods[method] += 1

    # Coverage and identity stats
    try:
        cov = float(row.get("% Coverage of reference", 0))
        coverages.append(cov)
    except (ValueError, TypeError):
        pass
    try:
        ident = float(row.get("% Identity to reference", 0))
        identities.append(ident)
    except (ValueError, TypeError):
        pass

total_genes = len(rows)
amr_genes = type_counts.get("AMR", 0)
stress_genes = type_counts.get("STRESS", 0)
virulence_genes = type_counts.get("VIRULENCE", 0)

# --- File 1: General Stats ---
with open(f"{sample}_amrfinder_generalstats_mqc.tsv", "w") as fh:
    fh.write("# id: 'amrfinder_generalstats'\\n")
    fh.write("# plot_type: 'generalstats'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     - Total_Genes:\\n")
    fh.write("#         title: 'AMRFinder Genes'\\n")
    fh.write("#         description: 'Total genes detected by AMRFinderPlus'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Blues'\\n")
    fh.write("#         placement: 2000\\n")
    fh.write("#     - AMR_Genes:\\n")
    fh.write("#         title: 'AMR Genes'\\n")
    fh.write("#         description: 'Antimicrobial resistance genes'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Reds'\\n")
    fh.write("#         placement: 2001\\n")
    fh.write("#     - Stress_Genes:\\n")
    fh.write("#         title: 'Stress Genes'\\n")
    fh.write("#         description: 'Stress response genes'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Greens'\\n")
    fh.write("#         placement: 2002\\n")
    fh.write("#     - Virulence_Genes:\\n")
    fh.write("#         title: 'Virulence Genes'\\n")
    fh.write("#         description: 'Virulence factor genes'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Purples'\\n")
    fh.write("#         placement: 2003\\n")
    fh.write("Sample\\tTotal_Genes\\tAMR_Genes\\tStress_Genes\\tVirulence_Genes\\n")
    fh.write(f"{sample}\\t{total_genes}\\t{amr_genes}\\t{stress_genes}\\t{virulence_genes}\\n")

# --- File 2: AMR Class Distribution (Bargraph) ---
with open(f"{sample}_amrfinder_classes_mqc.tsv", "w") as fh:
    fh.write("# id: 'amrfinder_amr_classes'\\n")
    fh.write("# section_name: 'AMRFinder: AMR Drug Classes'\\n")
    fh.write("# description: 'Distribution of antimicrobial resistance gene classes detected by AMRFinderPlus.'\\n")
    fh.write("# plot_type: 'bargraph'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     title: 'AMR Gene Classes'\\n")
    fh.write("#     ylab: 'Gene Count'\\n")
    fh.write("#     stacking: 'normal'\\n")

    if amr_classes:
        sorted_classes = sorted(amr_classes.keys())
        fh.write("Sample\\t" + "\\t".join(sorted_classes) + "\\n")
        counts = [str(amr_classes[c]) for c in sorted_classes]
        fh.write(f"{sample}\\t" + "\\t".join(counts) + "\\n")
    else:
        fh.write("Sample\\tNo_AMR_Detected\\n")
        fh.write(f"{sample}\\t0\\n")

# --- File 3: Detail Table ---
avg_cov = sum(coverages) / len(coverages) if coverages else 0.0
avg_ident = sum(identities) / len(identities) if identities else 0.0

# Standard AMRFinder methods to always show columns for
standard_methods = ["EXACTX", "BLASTX", "ALLELEX", "PARTIALX",
                    "PARTIAL_CONTIG_ENDX", "INTERNAL_STOP", "HMM"]

with open(f"{sample}_amrfinder_detail_mqc.tsv", "w") as fh:
    fh.write("# id: 'amrfinder_detail'\\n")
    fh.write("# section_name: 'AMRFinder: Detection Methods'\\n")
    fh.write("# description: 'AMRFinderPlus detection method breakdown and quality metrics per sample.'\\n")
    fh.write("# plot_type: 'table'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     id: 'amrfinder_detail_table'\\n")
    fh.write("#     title: 'AMRFinder Detection Details'\\n")
    fh.write("#     namespace: 'AMRFinder'\\n")

    # Build header: methods + quality metrics
    method_cols = [m for m in standard_methods if methods.get(m, 0) > 0]
    # Also include any non-standard methods found
    for m in sorted(methods.keys()):
        if m not in standard_methods and m not in method_cols:
            method_cols.append(m)

    # If no methods found (empty AMRFinder output), add placeholder
    if not method_cols:
        method_cols = ["No_Genes_Detected"]
        methods["No_Genes_Detected"] = 0

    header_parts = ["Sample"] + method_cols + ["Avg_Coverage_Pct", "Avg_Identity_Pct"]
    fh.write("\\t".join(header_parts) + "\\n")

    data_parts = [sample]
    for m in method_cols:
        data_parts.append(str(methods.get(m, 0)))
    data_parts.append(f"{avg_cov:.1f}")
    data_parts.append(f"{avg_ident:.1f}")
    fh.write("\\t".join(data_parts) + "\\n")

PYTHON_SCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
