/*
 * EXTRACT_CODON_FREQS - Per-read codon extraction and AA frequency analysis
 *
 * Uses pysam to extract 3-base codon triplets from individual reads at
 * user-specified protein positions, translates to amino acids, and computes
 * frequency distributions. This read-level approach correctly handles
 * codon linkage (all 3 bases from the same read) — unlike per-base mpileup
 * which would lose the codon structure.
 *
 * Designed for pooled/metagenomic samples where multiple alleles may coexist.
 * Reports frequency distributions rather than consensus calls.
 *
 * Input CDS FASTA assumption: position 1 = first codon (ATG start).
 * No UTRs, no introns.
 */

process EXTRACT_CODON_FREQS {

    label 'ultralow'
    container 'quay.io/biocontainers/mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45:3e0fc1ebdf2007459f18c33c65d38d2b031b0052-0'

    input:
        tuple val(sample), path(bam), path(bai)
        path(cds_fasta)
        path(positions_csv)

    output:
        tuple val(sample), path("${sample}_codon_freqs.tsv"), emit: codon_freqs
        path("versions.yml"),                                  emit: versions

    script:
    """
    python3 << 'PYTHON_SCRIPT'
import pysam
import csv
from collections import Counter

# Standard genetic code translation table
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

sample = "${sample}"
min_base_qual = int(${params.snp_min_base_quality})
min_depth = int(${params.snp_min_depth})

# -------------------------------------------------------
# 1. Read positions of interest from CSV
# -------------------------------------------------------
positions = []
with open("${positions_csv}") as f:
    reader = csv.DictReader(f)
    for row in reader:
        positions.append(row)

# -------------------------------------------------------
# 2. Open BAM file
# -------------------------------------------------------
bam = pysam.AlignmentFile("${bam}", "rb")
bam_references = set(bam.references)

# -------------------------------------------------------
# 3. For each position, extract codon frequencies from reads
# -------------------------------------------------------
results = []
for pos in positions:
    gene = pos['gene_name']
    prot_pos = int(pos['protein_position'])
    wt_aa = pos['wt_aa']
    mutant_aas = pos.get('mutant_aas', '')
    annotation = pos.get('annotation', '')

    # CDS coordinates (0-based):
    # Protein position P maps to nucleotides (P-1)*3 to (P-1)*3+2
    codon_start = (prot_pos - 1) * 3   # 0-based inclusive
    codon_end = codon_start + 3         # 0-based exclusive

    aa_counts = Counter()
    codon_counts = Counter()

    # If gene not in BAM references, no reads aligned — report depth=0
    if gene not in bam_references:
        results.append({
            'gene': gene, 'pos': prot_pos, 'wt': wt_aa,
            'depth': 0, 'aa_dist': '', 'codons': '',
            'mutants_expected': mutant_aas, 'mutants_found': '',
            'annotation': annotation
        })
        continue

    # Iterate reads overlapping this codon region
    for read in bam.fetch(gene, codon_start, codon_end):
        # Skip unmapped, secondary, and supplementary alignments
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Get aligned pairs: list of (read_pos, ref_pos)
        # read_pos is None for deletions; ref_pos is None for insertions
        aligned_pairs = read.get_aligned_pairs()
        codon_bases = {}
        codon_quals = {}

        for read_pos, ref_pos in aligned_pairs:
            if ref_pos is not None and codon_start <= ref_pos < codon_end:
                offset = ref_pos - codon_start  # 0, 1, or 2
                if read_pos is not None:
                    codon_bases[offset] = read.query_sequence[read_pos]
                    codon_quals[offset] = read.query_qualities[read_pos]

        # Only count reads where all 3 codon bases are present and pass quality
        if len(codon_bases) == 3 and all(
            codon_quals.get(i, 0) >= min_base_qual for i in range(3)
        ):
            codon = (codon_bases[0] + codon_bases[1] + codon_bases[2]).upper()
            aa = CODON_TABLE.get(codon, '?')
            aa_counts[aa] += 1
            codon_counts[codon] += 1

    # Format AA frequency distribution (descending by count)
    total = sum(aa_counts.values())
    aa_dist_parts = []
    for aa, count in sorted(aa_counts.items(), key=lambda x: -x[1]):
        freq = count / total if total > 0 else 0
        aa_dist_parts.append(f"{aa}:{freq:.4f}")

    # Format codon frequency distribution
    codon_dist_parts = []
    for codon, count in sorted(codon_counts.items(), key=lambda x: -x[1]):
        freq = count / total if total > 0 else 0
        codon_dist_parts.append(f"{codon}:{freq:.4f}")

    # Check which expected mutant AAs were actually observed
    expected = [m.strip() for m in mutant_aas.split(';') if m.strip()]
    found = [m for m in expected if m in aa_counts]

    results.append({
        'gene': gene, 'pos': prot_pos, 'wt': wt_aa,
        'depth': total, 'aa_dist': ','.join(aa_dist_parts),
        'codons': ','.join(codon_dist_parts),
        'mutants_expected': mutant_aas, 'mutants_found': ';'.join(found),
        'annotation': annotation
    })

bam.close()

# -------------------------------------------------------
# 4. Write per-position output TSV
# -------------------------------------------------------
with open(f"{sample}_codon_freqs.tsv", 'w') as f:
    header = [
        'sample', 'gene_name', 'protein_position', 'wt_aa', 'codon_depth',
        'aa_distribution', 'codon_distribution', 'mutant_aas_expected',
        'mutants_detected', 'annotation'
    ]
    f.write('\\t'.join(header) + '\\n')
    for r in results:
        row = [
            sample, r['gene'], str(r['pos']), r['wt'], str(r['depth']),
            r['aa_dist'], r['codons'], r['mutants_expected'],
            r['mutants_found'], r['annotation']
        ]
        f.write('\\t'.join(row) + '\\n')

PYTHON_SCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
