#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p long
#SBATCH -o "balrog.stdout.%j.%N"
#SBATCH -e "balrog.stderr.%j.%N"
#SBATCH --mail-user=your.email@usda.gov
#SBATCH --mail-type=START,END,FAIL

module load apptainer
# Requires Nextflow >=24.04 (process.resourceLimits directive).
# Run `module avail nextflow` and pick the newest >=24.04 build available on Ceres.
module load nextflow/24.10.4

nextflow run main.nf \
    -c configs/ceres/ceres.cfg \
    --slurm_account your_group \
    --sample_sheet samplesheet.csv \
    --kraken2_db /path/to/kraken2_db \
    --sylph_db /path/to/sylph.syldb \
    --diamond_db /path/to/amr.dmnd \
    --outdir ./results \
    -resume
