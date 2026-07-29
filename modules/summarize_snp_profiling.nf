/*
 * SUMMARIZE_SNP_PROFILING
 *
 * Reads per-sample codon frequency TSV and produces MultiQC Custom Content files:
 *   1. General stats  (_snp_generalstats_mqc.tsv)  — summary row per sample
 *   2. Detail table   (_snp_detail_mqc.tsv)        — one row per gene/position
 *
 * These _mqc.tsv files are auto-detected by MultiQC's Custom Content module.
 * Follows the same pattern as summarize_amrfinder.nf.
 */

process SUMMARIZE_SNP_PROFILING {

    label 'ultralow'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    input:
        tuple val(sample), path(codon_freqs_tsv)

    output:
        path("${sample}_snp_generalstats_mqc.tsv"), emit: generalstats
        path("${sample}_snp_detail_mqc.tsv"),       emit: detail
        path("versions.yml"),                        emit: versions

    script:
    """
    python3 << 'PYTHON_SCRIPT'
import csv

sample = "${sample}"
tsv_path = "${codon_freqs_tsv}"
min_depth = int(${params.snp_min_depth})

# --- Read codon frequency data ---
rows = []
with open(tsv_path, "r") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        rows.append(row)

# --- Compute summary statistics ---
total_positions = len(rows)
mutants_detected = 0
depths = []

for row in rows:
    depth = int(row.get('codon_depth', 0))
    depths.append(depth)
    detected = row.get('mutants_detected', '').strip()
    if detected:
        mutants_detected += 1

avg_depth = sum(depths) / len(depths) if depths else 0.0
min_depth_val = min(depths) if depths else 0

# --- File 1: General Stats ---
with open(f"{sample}_snp_generalstats_mqc.tsv", "w") as fh:
    fh.write("# id: 'snp_profiling_generalstats'\\n")
    fh.write("# plot_type: 'generalstats'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     - Positions_Profiled:\\n")
    fh.write("#         title: 'SNP Positions'\\n")
    fh.write("#         description: 'Number of target positions profiled'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Blues'\\n")
    fh.write("#         placement: 3000\\n")
    fh.write("#     - Mutants_Detected:\\n")
    fh.write("#         title: 'Mutants Found'\\n")
    fh.write("#         description: 'Positions where known mutant AAs were observed'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Reds'\\n")
    fh.write("#         placement: 3001\\n")
    fh.write("#     - Avg_Depth:\\n")
    fh.write("#         title: 'SNP Avg Depth'\\n")
    fh.write("#         description: 'Average codon read depth across profiled positions'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Greens'\\n")
    fh.write("#         format: '{:,.1f}'\\n")
    fh.write("#         placement: 3002\\n")
    fh.write("#     - Min_Depth:\\n")
    fh.write("#         title: 'SNP Min Depth'\\n")
    fh.write("#         description: 'Minimum codon read depth across profiled positions'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Oranges'\\n")
    fh.write("#         placement: 3003\\n")
    fh.write("Sample\\tPositions_Profiled\\tMutants_Detected\\tAvg_Depth\\tMin_Depth\\n")
    fh.write(f"{sample}\\t{total_positions}\\t{mutants_detected}\\t{avg_depth:.1f}\\t{min_depth_val}\\n")

# --- File 2: Detail Table ---
with open(f"{sample}_snp_detail_mqc.tsv", "w") as fh:
    fh.write("# id: 'snp_profiling_detail'\\n")
    fh.write("# section_name: 'SNP Profiling: Variant Details'\\n")
    fh.write("# description: 'Per-position amino acid frequency distributions at targeted SNP sites. "
             "Each row shows the observed AA distribution at a specific codon position.'\\n")
    fh.write("# plot_type: 'table'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     id: 'snp_profiling_detail_table'\\n")
    fh.write("#     title: 'Targeted SNP Profiling'\\n")
    fh.write("#     namespace: 'SNP Profiling'\\n")
    fh.write("#     no_beeswarm: true\\n")

    # Header
    fh.write("Sample\\tGene\\tPosition\\tWT_AA\\tDepth\\tAA_Distribution\\t"
             "Mutants_Expected\\tMutants_Detected\\tAnnotation\\n")

    for row in rows:
        gene = row.get('gene_name', '')
        pos = row.get('protein_position', '')
        wt = row.get('wt_aa', '')
        depth = row.get('codon_depth', '0')
        aa_dist = row.get('aa_distribution', '')
        expected = row.get('mutant_aas_expected', '')
        detected = row.get('mutants_detected', '')
        annotation = row.get('annotation', '')

        # Add low-depth warning to annotation
        if int(depth) < min_depth and int(depth) > 0:
            annotation = f"[LOW DEPTH] {annotation}" if annotation else "[LOW DEPTH]"
        elif int(depth) == 0:
            annotation = f"[NO COVERAGE] {annotation}" if annotation else "[NO COVERAGE]"

        # Row ID: {sample}__{gene}_{position} for unique rows in merged table
        row_id = f"{sample}__{gene}_{pos}"

        # Format detected: show 'none' if expected mutants exist but none found
        if expected and not detected:
            detected = "none"

        fh.write(f"{row_id}\\t{gene}\\t{pos}\\t{wt}\\t{depth}\\t{aa_dist}\\t"
                 f"{expected}\\t{detected}\\t{annotation}\\n")

PYTHON_SCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
