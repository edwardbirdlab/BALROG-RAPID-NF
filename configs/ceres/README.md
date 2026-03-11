### Running BALROG-RAPID on USDA ARS Ceres

> This file is meant to help you successfully run the BALROG-RAPID pipeline on Ceres SCINet services.

```bash
# From the BALROG-RAPID root directory
cp nextflow/configs/ceres/balrog.template.slurm balrog.slurm
vim balrog.slurm      # edit with your email and database paths
vim samplesheet.csv   # create samplesheet (columns: sample,r1,r2)
sbatch balrog.slurm
```

Optional: to run host profiling, add `--host_sheet hosts.csv` to the nextflow command in the slurm script.

This file and all files in this directory are in the public domain.
