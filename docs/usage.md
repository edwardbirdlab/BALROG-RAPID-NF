# Usage

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 24.04 (the pipeline uses the `process.resourceLimits` directive)
- One of Docker, Singularity, or Apptainer

## Container engine

No container engine is enabled by default. Select one explicitly with `-profile`:

```bash
nextflow run main.nf -profile docker    ...
nextflow run main.nf -profile singularity ...
nextflow run main.nf -profile apptainer ...
```

Running without a `-profile` flag executes every tool directly on `$PATH`, which will fail unless you have
all the underlying tools installed locally -- always pass a profile.

## Samplesheet

A single unified CSV drives both short-read (Illumina/Element) and long-read (Oxford Nanopore) samples:

```csv
sample,library,molecule,r1,r2,reads
Sample1,Short,gDNA,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz,
Sample2,Long,gDNA,,,/path/to/Sample2.fastq.gz
```

- `library`: `Short` (uses `r1`/`r2`) or `Long` (uses `reads`)
- `molecule`: free-text, carried through for your own record-keeping
- Rows sharing a `sample` name are concatenated (re-sequencing support) as long as `library`/`molecule`
  are consistent across those rows
- Short- and long-read samples can be mixed in one samplesheet; each routes to the matching workflow and
  both merge into a single MultiQC report

## Required parameters

| Parameter | Description |
|-----------|--------------|
| `--sample_sheet` | Unified CSV: `sample,library,molecule,r1,r2,reads` |
| `--kraken2_db` | Path to a unified Kraken2 database (required when `--run_taxonomy`, on by default) |
| `--sylph_db` | Path to a Sylph database (`.syldb`) (required when `--run_taxonomy`) |
| `--diamond_db` | Path to a Diamond AMR database (`.dmnd`) (required when `--run_amr`, on by default) |

## Optional steps

| Flag | Adds | Notes |
|------|------|-------|
| `--host_sheet hosts.csv` | Kraken2 host-contamination profiling | CSV: `host_name,kraken2_db`. Every sample is classified against every host DB; reads are not filtered. |
| `--run_bbduk` | BBDuk adapter trimming (Element Aviti, short reads) | Runs after FASTP, before the trimmed FastQC. `--bbduk_adapters` defaults to the bundled `assets/element_aviti_adapters.fasta`. |
| `--run_spike_in` | Bowtie2/samtools T. thermophilus spike-in removal (short reads only) | Runs after QC, before taxonomy/host/AMR. `--spike_in_bt2` defaults to the bundled `assets/t_thermophilus_bt2/` index. |
| `--run_nonpareil` (on by default) | Nonpareil coverage estimation | Runs on Kraken2-extracted bacterial reads when `--run_taxonomy` is also on; otherwise on the full read set. Disable with `--run_nonpareil false`. |
| `--custom_qc` | `% Bacterial reads` column in the MultiQC General Stats table | Requires `--run_taxonomy`. |
| `--run_snp_profiling` | Targeted amino-acid variant profiling at specific CDS positions | Requires `--snp_cds_fasta` and `--snp_positions_csv`. BWA-MEM (short reads) or minimap2 (long reads) alignment, per-read codon extraction via pysam. |

Long-read samples automatically use Porechop_ABI + Filtlong + NanoPlot for QC (instead of FASTP/BBDuk) and
Flye for AMR-subset assembly (instead of SPAdes) -- no extra flags needed, this follows from `library=Long`
in the samplesheet. `--filtlong_min_length`, `--filtlong_keep_percent`, and `--flye_read_type` tune these steps.

See the top-level [README](../README.md) for full parameter tables, a pipeline diagram, and per-feature
usage examples.

## Ceres (USDA ARS SCINet)

See [`configs/ceres/README.md`](../configs/ceres/README.md) and the
[`configs/ceres/balrog.template.slurm`](../configs/ceres/balrog.template.slurm) launcher template.
`--slurm_account` is required when using `-c configs/ceres/ceres.cfg`.
