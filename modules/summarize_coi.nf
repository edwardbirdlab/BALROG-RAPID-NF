/*
 * SUMMARIZE_COI
 *
 * "The pipeline is: run KMA, filter two columns." Filters KMA's .res to
 * calls where BOTH Template_Coverage and Query_Identity clear their floors
 * (neither works alone on COI -- coverage alone lets a wrong-family
 * reference through at 100% coverage/84% identity because ~3 short reads
 * tile a 658bp barcode; identity alone lets a single conserved-region read
 * through at 100% identity/3% coverage). No tiers, no rank-by-identity --
 * below species, rank is not recoverable from percent identity on COI (see
 * docs/usage.md).
 *
 * Parses the Template field positionally (column order, not names --
 * KMA/BOLD header styles are not stable across versions):
 *   "<taxid>|<BIN>_<processid> <taxon>"   e.g. "10004701|BOLD:AAN5248_PHNXG853-13 Cecidomyiidae"
 *   "<taxid>|<processid> <taxon>"         when the reference has no BIN
 *
 * A sample with zero calls is a real result (this is common -- roughly 1 in
 * 6-12 real libraries in validation), not an error: the detail table still
 * gets a "no confident hits" row rather than being empty or skipped.
 *
 * Optionally joins a BIN/processid -> full lineage lookup (built once by
 * MAKE_LINEAGE_LOOKUP) when --coi_lineage_table is set; unmatched or
 * unavailable lineage is reported as "?", never silently blank.
 *
 * Emits MultiQC Custom Content files:
 *   1. General stats  (_coi_generalstats_mqc.tsv)
 *   2. Detail table    (_coi_detail_mqc.tsv)
 */

process SUMMARIZE_COI {

    label 'ultralow'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    input:
        tuple val(sample), path(res)
        path(lineage_lookup)

    output:
        path("${sample}_coi_generalstats_mqc.tsv"), emit: generalstats
        path("${sample}_coi_detail_mqc.tsv"),        emit: detail
        path("versions.yml"),                         emit: versions

    script:
    """
    python3 << 'PYTHON_SCRIPT'
sample = "${sample}"
min_coverage = float(${params.coi_min_coverage})
min_identity = float(${params.coi_min_identity})


def parse_template(template):
    # BIN/processid/taxon, split by position -- see module header.
    taxid, sep, rest = template.partition("|")
    token, sep, taxon = rest.partition(" ")
    bin_id, pid = "", token
    if token.startswith("BOLD:"):
        bin_id, sep, pid = token.partition("_")
    return bin_id, pid, taxon


# --- Optional lineage lookup (processid -> lineage string) ---
lineage_by_pid = {}
with open("${lineage_lookup}") as fh:
    header = fh.readline()
    for line in fh:
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) >= 3:
            lineage_by_pid[parts[0]] = parts[2]

# --- Filter KMA's .res: BOTH coverage and identity must clear their floors ---
calls = []
with open("${res}") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 9:
            continue
        coverage = float(cols[5])
        identity = float(cols[6])
        if coverage < min_coverage or identity < min_identity:
            continue
        depth = cols[8].strip()
        bin_id, pid, taxon = parse_template(cols[0].strip())
        calls.append({
            "bin": bin_id,
            "processid": pid,
            "taxon": taxon.strip(),
            "coverage": coverage,
            "identity": identity,
            "depth": depth,
            "lineage": lineage_by_pid.get(pid, "?"),
        })

calls.sort(key=lambda c: c["coverage"], reverse=True)

n_calls = len(calls)
top_taxon = "none"
if calls:
    top = calls[0]
    top_taxon = top["taxon"] or top["bin"] or top["processid"]

# --- File 1: General Stats ---
with open(f"{sample}_coi_generalstats_mqc.tsv", "w") as fh:
    fh.write("# id: 'coi_generalstats'\\n")
    fh.write("# plot_type: 'generalstats'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     - COI_N_Calls:\\n")
    fh.write("#         title: 'COI Calls'\\n")
    fh.write("#         description: 'Reference calls clearing both the coverage and identity floors (COI/KMA)'\\n")
    fh.write("#         min: 0\\n")
    fh.write("#         scale: 'Greens'\\n")
    fh.write("#         placement: 4000\\n")
    fh.write("#     - COI_Top_Taxon:\\n")
    fh.write("#         title: 'Top COI Call'\\n")
    fh.write("#         description: 'Highest-coverage COI call (taxon name, or BIN/processid when unnamed)'\\n")
    fh.write("#         placement: 4001\\n")
    fh.write("Sample\\tCOI_N_Calls\\tCOI_Top_Taxon\\n")
    fh.write(f"{sample}\\t{n_calls}\\t{top_taxon}\\n")

# --- File 2: Detail Table ---
with open(f"{sample}_coi_detail_mqc.tsv", "w") as fh:
    fh.write("# id: 'coi_detail'\\n")
    fh.write("# section_name: 'COI Insect ID'\\n")
    fh.write("# description: 'KMA calls against the BOLD COI Diptera reference set, passing both the "
             "coverage and identity floors. BIN is reported alongside/instead of species name -- most "
             "BOLD COI records carry no species-level name. A sample with no rows below still passed "
             "through this step; it means no reference cleared both floors, not that the step failed.'\\n")
    fh.write("# plot_type: 'table'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     id: 'coi_detail_table'\\n")
    fh.write("#     title: 'COI Insect ID'\\n")
    fh.write("#     namespace: 'COI Insect ID'\\n")

    fh.write("Sample\\tBIN\\tProcessID\\tTaxon\\tCoverage\\tIdentity\\tDepth\\tLineage\\n")

    if not calls:
        fh.write(f"{sample}\\t-\\t-\\tno confident hits\\t0\\t0\\t0\\t?\\n")
    else:
        for i, c in enumerate(calls):
            row_id = f"{sample}__{i}__{c['processid']}"
            fh.write(f"{row_id}\\t{c['bin']}\\t{c['processid']}\\t{c['taxon']}\\t"
                     f"{c['coverage']:.2f}\\t{c['identity']:.2f}\\t{c['depth']}\\t{c['lineage']}\\n")

PYTHON_SCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
