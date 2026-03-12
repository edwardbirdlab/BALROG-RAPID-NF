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
  ├─ 1. Read QC (FASTP + FastQC)
  │
  ├─ 2. Taxonomy (Kraken2 + Sylph)        ─┐
  ├─ 3. Host Profiling (Kraken2 x N hosts)  ├─ run in parallel
  └─ 4. AMR Detection                      ─┘
         Diamond BLASTX → seqtk subseq → SPAdes micro-assembly → AMRFinder
  │
  ├─ 5. Software Versions (collected from all steps)
  └─ 6. MultiQC (aggregated QC report)
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
│   └── fastqc_trimmed/
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
| MultiQC | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |
