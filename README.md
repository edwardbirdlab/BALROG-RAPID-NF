# BALROG-RAPID Nextflow Pipeline

Nextflow DSL2 pipeline for metagenomic pathogen detection, host profiling, and AMR detection. This pipeline was developed with the intention of producing rapid results, formatted for downstream visualization with the BALROG-RAPID-DASHBOARD package (link).

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

Add pipeline flowchart here

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

Each sample is classified against every host database. This is for contamination reporting only; meaning reads are **not** filtered. This is meant as a QC step to see host contamination.

## Element Biosciences Aviti Adapter Trimming (Optional)

When sequencing on an Element Biosciences Aviti instrument, reads may contain platform-specific sequencing artifacts. Enable the optional BBDuk trimming step to clean these:

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --run_bbduk \
  ...
```

BBDuk runs **after FASTP, but before FastQC (Trimmed)**, so the final QC reflects fully cleaned reads. The adapter FASTA defaults to the Element Aviti concatenated adapter file and does not need to be specified unless you have a custom one.

The adapter FASTA is bundled in the pipeline at `assets/element_aviti_adapters.fasta` .

## T. thermophilus Spike-in Removal (Optional)

When T. thermophilus is used as a spike-in control, enable Bowtie2/Samtools to remove spike-in reads before downstream analysis:

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --run_spike_in \
  ...
```

A pre-built Bowtie2 index for T. thermophilus is bundled in the pipeline at `assets/t_thermophilus_bt2/`. The step runs **after QC but before taxonomy, host profiling, and AMR detection**, so all downstream analyses use spike-depleted reads.

To use a custom spike-in reference, build a Bowtie2 index and pass the directory:

```bash
--spike_in_bt2 /path/to/custom_bt2_index/
```

## Nonpareil Coverage Estimation (Optional)

Nonpareil estimates metagenomic sequencing coverage by analyzing read redundancy.

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --run_nonpareil \
  ...
```

**Bacterial read filtering**: When taxonomy profiling is also enabled (`--run_taxonomy`, default true), Nonpareil automatically runs on **bacterial reads only**. Kraken2 output is used to extract reads classified under Bacteria (taxid 2, including all children) via KrakenTools. This removes the host contamination signal that would otherwise bias coverage estimates. Extraction statistics are published to `results/nonpareil/bacterial_extraction/`.

When taxonomy is disabled, Nonpareil runs on the full read set (not recommended unless there is very little host contamination).

Nonpareil runs on R1 reads only (one read per pair) using the kmer algorithm. Coverage estimates appear in the MultiQC report and raw output files are published to `results/nonpareil/`.

## Custom QC Metrics (Optional)

Adds custom metrics to the MultiQC General Stats table. Currently reports the **percentage of reads classified as Bacteria** by Kraken2, may include normalized T. thermophilus normalized results in the future.

```bash
nextflow run nextflow/main.nf \
  --sample_sheet samplesheet.csv \
  --custom_qc \
  ...
```

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
| `--custom_qc` | false | Enable custom QC metrics (% bacterial reads) |
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
│   ├── {sample}.npa
│   └── bacterial_extraction/                      # only when taxonomy + nonpareil both enabled
│       └── {sample}_bacterial_extraction_stats.tsv
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
| Bowtie2 + Samtools (Spike-in) | `quay.io/biocontainers/mulled-v2-229691629e0b12c862d76101f90a597d5c1c81d4:484c804e1d5952c9023891b6f9a19f7f15815145-0` (Bowtie2 2.4.5 + Samtools 1.16.1) |
| KrakenTools | `quay.io/biocontainers/krakentools:1.2.1--pyh7e72e81_0` |
| Nonpareil | `quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2` |
| BBMap (BBDuk) | `ebird013/bbmap:latest` |
| Sylph-tax | `quay.io/biocontainers/sylph-tax:1.8.0--pyhdfd78af_0` |
| MultiQC | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |
