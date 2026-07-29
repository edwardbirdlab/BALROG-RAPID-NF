# Output

This document describes the output produced under `--outdir` (default `./results`). Short-read and
long-read samples publish into the same directory layout, keyed by `{sample}` rather than read type.

```
results/
├── qc/
│   ├── fastqc_raw/                                 # FastQC on raw reads (both read types)
│   ├── fastp/                                      # FASTP trimming report + trimmed reads (short reads)
│   ├── bbduk/                                      # short reads, only when --run_bbduk
│   ├── porechop_abi/                               # long reads: adapter trimming log + trimmed reads
│   ├── filtlong/                                   # long reads: length/quality filtering log + filtered reads
│   ├── nanoplot/                                   # long reads: QC plots + NanoStats.txt
│   └── fastqc_trimmed/                             # FastQC on final trimmed reads (both read types)
├── spike_in/                                       # short reads only, when --run_spike_in
│   ├── {sample}_spike_stats.tsv
│   └── {sample}_bowtie2_spike.log
├── nonpareil/                                      # only when --run_nonpareil (on by default)
│   ├── {sample}.npo
│   ├── {sample}.npa
│   └── bacterial_extraction/                       # only when taxonomy + nonpareil both enabled
│       └── {sample}_bacterial_extraction_stats.tsv
├── taxonomy/
│   ├── kraken2/{sample}_k2report.tsv                # unified-DB Kraken2 report
│   └── sylph/{sample}_profile.tsv                   # raw Sylph profile
├── host/
│   └── kraken2/{sample}_{host}_k2report.tsv         # one report per sample x host DB
├── snp_profiling/                                   # only when --run_snp_profiling
│   ├── alignments/{sample}_snp.bam                  # BWA-MEM (short) or minimap2 (long) alignment to CDS reference
│   └── codon_freqs/{sample}_codon_freqs.tsv          # per-position amino acid / codon frequency distributions
├── amr/
│   ├── diamond/{sample}_diamond_hits.tsv            # Diamond BLASTX hits (AMR candidate reads)
│   ├── extracted_reads/{sample}_amr_R{1,2}.fastq.gz # seqtk-subset reads (short reads)
│   ├── extracted_reads/{sample}_amr_reads.fastq.gz  # seqtk-subset reads (long reads)
│   ├── assemblies/{sample}_contigs.fasta            # SPAdes (short reads) or Flye (long reads) micro-assembly
│   └── amrfinder/{sample}_amrfinder.tsv             # AMRFinderPlus results
├── multiqc/
│   ├── multiqc_report.html                          # aggregated QC report (both read types merged)
│   └── multiqc_report_data/
└── pipeline_info/
    └── software_versions.yml                        # deduplicated tool versions from every process
```

## Key deliverables

- **`multiqc/multiqc_report.html`** -- the single-page QC summary: read trimming/filtering, taxonomy, host
  contamination, AMR gene detection, targeted SNP profiling, and (when enabled) coverage estimation and
  custom QC metrics, merged across short-read and long-read samples.
- **`taxonomy/kraken2/*_k2report.tsv`** and **`amr/amrfinder/*_amrfinder.tsv`** are the pipeline's primary
  per-sample taxonomy and AMR-gene deliverables.
- **`snp_profiling/codon_freqs/*_codon_freqs.tsv`** reports amino-acid frequency distributions at
  user-specified codon positions, for surveillance of known resistance/adaptive mutations in pooled samples.
- **`pipeline_info/software_versions.yml`** records the exact tool version used in every process, for
  a specific run's provenance.

## Notes

- Directories for optional steps (`spike_in/`, `nonpareil/`, `bbduk/`, `snp_profiling/`) only appear when
  the corresponding `--run_*` flag is set.
- Host profiling produces one Kraken2 report per `(sample, host)` pair when `--host_sheet` is provided;
  it never filters reads, so downstream analyses always run on host-inclusive read sets.
- Long-read samples never produce `spike_in/` or `fastp/`/`bbduk/` output -- those steps are short-read
  specific (see [usage.md](usage.md)).
