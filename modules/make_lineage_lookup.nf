/*
 * MAKE_LINEAGE_LOOKUP
 *
 * Builds a processid -> BOLD lineage lookup once (not per-sample) from the
 * corrected COI dereplication table, for SUMMARIZE_COI to join onto calls.
 * Optional: only runs when --coi_lineage_table is set.
 *
 * BOLD writes the literal string "None" for missing values -- truthy in
 * most languages, and has silently corrupted analyses in this project
 * before. Strip it (and "NA"/"null") explicitly rather than trusting
 * truthiness.
 *
 * Use the CORRECTED derep table. An uncorrected one has been observed to
 * carry a mislabelled reference (a fly filed as a grasshopper), which
 * produces a report row where the taxon name and the joined lineage
 * flatly contradict each other.
 */

process MAKE_LINEAGE_LOOKUP {

    label 'ultralow'
    container 'quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0'

    input:
        path(derep_table)

    output:
        path("bold_lineage.tsv"), emit: lookup
        path("versions.yml"),     emit: versions

    script:
    """
    python3 << 'PYTHON_SCRIPT'
import csv

MISSING = {"", "None", "NA", "null"}
RANKS = [("kingdom", "k"), ("phylum", "p"), ("class", "c"), ("order", "o"),
         ("family", "f"), ("genus", "g"), ("species", "s")]

with open("${derep_table}", newline="") as fh, \\
     open("bold_lineage.tsv", "w") as out:
    reader = csv.DictReader(fh, delimiter="\\t")
    out.write("processid\\tbin\\tlineage\\n")
    for row in reader:
        parts = []
        for column, prefix in RANKS:
            value = (row.get(column) or "").strip()
            if value not in MISSING:
                parts.append(f"{prefix}__{value}")
        lineage = ";".join(parts)
        bin_uri = (row.get("bin_uri") or "").strip()
        if bin_uri in MISSING:
            bin_uri = ""
        out.write(f"{row.get('processid','')}\\t{bin_uri}\\t{lineage}\\n")

PYTHON_SCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
