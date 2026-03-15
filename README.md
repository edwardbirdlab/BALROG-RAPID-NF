# BALROG-RAPID Nextflow Pipeline

Nextflow DSL2 pipeline for metagenomic pathogen detection, host profiling, and AMR detection. This is Phase 1 of the BALROG-RAPID system — outputs feed into the Phase 2 risk map dashboard.

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 22.10
- Docker or Singularity (all tools run in containers)

## Quick Start

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --kraken2_db /path/to/kraken2_db \
  --sylph_db /path/to/sylph.syldb \
  --diamond_db /path/to/amr.dmnd \
  --outdir ./results
```

## Samplesheet

The pipeline uses a simple CSV with three columns:

```csv
sample,r1,r2
Sample1,/path/to/Sample1_R1_001.fastq.gz,/path/to/Sample1_R2_001.fastq.gz
Sample2,/path/to/Sample2_R1_001.fastq.gz,/path/to/Sample2_R2_001.fastq.gz
```

A separate metadata CSV (sample, site, coordinates, etc.) is used by Phase 2 and is not needed here.

## Pipeline Steps

```
Raw FASTQ
  │
  ├─ 1. Read QC
  │     FASTP → [optional BBDuk] → FastQC
  │
  ├─ 1b. [optional] Spike-in Removal (Bowtie2 vs T. thermophilus)
  │
  ├─ 2. Taxonomy (Kraken2 + Sylph)           ─┐
  ├─ 3. Host Profiling (Kraken2 x N hosts)    │
  ├─ 4. AMR Detection                         ├─ run in parallel
  │      Diamond BLASTX → seqtk → SPAdes → AMRFinder │
  └─ 5. [optional] Nonpareil Coverage        ─┘
  │
  ├─ 6. Software Versions (collected from all steps)
  └─ 7. MultiQC (aggregated QC report)
```

## Host Profiling (Optional)

To run Kraken2 against one or more host genomes, provide a host sheet CSV:

```bash
--host_sheet hosts.csv
```

```csv
host_name,kraken2_db
human,/path/to/kraken2_human_db
chicken,/path/to/kraken2_chicken_db
```

Each sample is classified against every host database. This is for contamination reporting only — reads are **not** filtered.

## Element Biosciences Aviti Adapter Trimming (Optional)

When sequencing on an Element Biosciences Aviti instrument, reads may contain platform-specific adapter sequences that FASTP doesn't fully remove. Enable the optional BBDuk trimming step to clean these:

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --run_bbduk \
  ...
```

BBDuk runs **after FASTP but before FastQC (Trimmed)**, so the final QC reflects fully cleaned reads. The adapter FASTA defaults to the Element Aviti concatenated adapter file and does not need to be specified unless you have a custom one.

BBDuk trimming stats appear automatically in the MultiQC report when enabled.

The adapter FASTA is bundled in the pipeline at `assets/element_aviti_adapters.fasta` and does not require internet access at runtime.

## T. thermophilus Spike-in Removal (Optional)

When T. thermophilus is used as a spike-in control, enable Bowtie2 alignment to remove spike-in reads before downstream analysis:

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --run_spike_in \
  ...
```

A pre-built Bowtie2 index for T. thermophilus is bundled in the pipeline at `assets/t_thermophilus_bt2/`. The step runs **after QC but before taxonomy, host profiling, and AMR detection**, so all downstream analyses use spike-depleted reads.

Reads are aligned with Bowtie2 and piped directly through samtools to extract unmapped pairs — any read pair where **either mate** maps to the spike-in genome is removed. This is more thorough than concordant-only filtering and avoids writing large intermediate SAM files to disk.

Spike-in alignment statistics appear in the MultiQC report and per-sample stats TSV files are published to `results/spike_in/`.

To use a custom spike-in reference, build a Bowtie2 index and pass the directory:

```bash
--spike_in_bt2 /path/to/custom_bt2_index/
```

## Nonpareil Coverage Estimation (Optional)

Nonpareil estimates metagenomic sequencing coverage by analyzing read redundancy — it answers "have we sequenced enough?" for each sample.

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --run_nonpareil \
  ...
```

Nonpareil runs on R1 reads only (one read per pair) using the kmer algorithm. It runs **in parallel** with taxonomy, host profiling, and AMR detection. Coverage estimates appear in the MultiQC report and raw output files are published to `results/nonpareil/`.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sample_sheet` | required | CSV with sample, r1, r2 columns |
| `--kraken2_db` | required | Path to unified Kraken2 database |
| `--sylph_db` | required | Path to Sylph database (.syldb) |
| `--diamond_db` | required | Path to Diamond AMR database (.dmnd) |
| `--outdir` | `./results` | Output directory |
| `--host_sheet` | null | Optional CSV of host databases |
| `--k2_confidence` | 0.3 | Kraken2 confidence threshold |
| `--k2_min_hit_groups` | 3 | Kraken2 minimum hit groups |
| `--fastp_q` | 20 | FASTP minimum quality score |
| `--fastp_minlen` | 100 | FASTP minimum read length |
| `--diamond_evalue` | 1e-10 | Diamond E-value threshold |
| `--diamond_max_targets` | 500 | Diamond max target sequences |
| `--run_qc` | true | Enable/disable read QC |
| `--run_bbduk` | false | Enable BBDuk adapter trimming (Element Aviti) |
| `--bbduk_adapters` | bundled | Path to adapter FASTA for BBDuk |
| `--bbduk_ktrim` | `r` | BBDuk kmer trim direction |
| `--bbduk_hdist` | 1 | BBDuk Hamming distance tolerance |
| `--bbduk_additional_args` | `''` | Extra bbduk.sh arguments |
| `--run_spike_in` | false | Enable T. thermophilus spike-in removal |
| `--spike_in_bt2` | bundled | Path to pre-built Bowtie2 index directory |
| `--run_nonpareil` | false | Enable Nonpareil coverage estimation |
| `--run_taxonomy` | true | Enable/disable taxonomy profiling |
| `--run_host_profiling` | true | Enable/disable host profiling |
| `--run_amr` | true | Enable/disable AMR detection |
| `--run_multiqc` | true | Enable/disable MultiQC report |
| `--multiqc_config` | built-in | Path to custom MultiQC config YAML |
| `--max_cpus` | 36 | Maximum CPUs per process |
| `--max_memory` | 128.GB | Maximum memory per process |

## Output Structure

```
results/
├── qc/
│   ├── fastqc_raw/
│   ├── fastp/
│   ├── bbduk/                                   # only when --run_bbduk
│   └── fastqc_trimmed/
├── spike_in/                                      # only when --run_spike_in
│   ├── {sample}_spike_stats.tsv
│   └── {sample}_bowtie2_spike.log
├── nonpareil/                                     # only when --run_nonpareil
│   ├── {sample}.npo
│   └── {sample}.npa
├── taxonomy/
│   ├── kraken2/{sample}_k2report.tsv
│   └── sylph/{sample}_profile.tsv
├── host/
│   └── kraken2/{sample}_{host}_k2report.tsv
├── amr/
│   ├── diamond/{sample}_diamond_hits.tsv
│   ├── extracted_reads/{sample}_amr_R{1,2}.fastq.gz
│   ├── assemblies/{sample}_contigs.fasta
│   └── amrfinder/{sample}_amrfinder.tsv
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
└── pipeline_info/
    └── software_versions.yml
```

The `taxonomy/` and `amr/amrfinder/` outputs are directly compatible with Phase 2 (`tools/balrog_riskmap.py`).

## Containers

All processes run in containers. No local tool installation needed.

| Tool | Container |
|------|-----------|
| FastQC | `ebird013/fastqc:0.12.1_custom` |
| FASTP | `biocontainers/fastp:v0.20.1_cv1` |
| Kraken2 | `ebird013/kraken2:2.1.3` |
| Sylph | `quay.io/biocontainers/sylph:0.9.0--ha6fb395_0` |
| Diamond | `quay.io/biocontainers/diamond:2.1.24--hf93d47f_0` |
| seqtk | `quay.io/biocontainers/seqtk:1.5--h577a1d6_1` |
| SPAdes | `ebird013/spades:3.15.5` |
| AMRFinder | `ncbi/amr:latest` |
| Bowtie2 + Samtools (Spike-in) | `quay.io/biocontainers/mulled-v2-...` (Bowtie2 2.4.5 + Samtools 1.16.1) |
| Nonpareil | `quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2` |
| BBMap (BBDuk) | `ebird013/bbmap:latest` |
| Sylph-tax | `quay.io/biocontainers/sylph-tax:1.8.0--pyhdfd78af_0` |
| MultiQC | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |
