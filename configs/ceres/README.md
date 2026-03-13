### Running BALROG-RAPID on USDA ARS Ceres

> This file is meant to help you successfully run the BALROG-RAPID pipeline on Ceres SCINet services.

```bash
# From the BALROG-RAPID root directory
cp nextflow/configs/ceres/balrog.template.slurm balrog.slurm
vim balrog.slurm      # edit with your email, database paths, and SLURM account
vim samplesheet.csv   # create samplesheet (columns: sample,r1,r2)
sbatch balrog.slurm
```

#### SLURM Account (`--slurm_account`)

The `--slurm_account` parameter is **required** when using the Ceres config. It specifies the SLURM allocation group for job billing (the `-A` flag). The pipeline will error at job submission if this is not set.

```bash
nextflow run nextflow/main.nf \
    -c nextflow/configs/ceres/ceres.cfg \
    --slurm_account your_group \
    --sample_sheet samplesheet.csv \
    ...
```

Optional: to run host profiling, add `--host_sheet hosts.csv` to the nextflow command in the slurm script.

This file and all files in this directory are in the public domain.
