#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p long
#SBATCH -o "balrog.stdout.%j.%N"
#SBATCH -e "balrog.stderr.%j.%N"
#SBATCH --mail-user=your.email@usda.gov
#SBATCH --mail-type=START,END,FAIL

module load apptainer
module load nextflow/23.10.1

nextflow run nextflow/main.nf \
    -c nextflow/configs/ceres/ceres.cfg \
    --sample_sheet samplesheet.csv \
    --kraken2_db /path/to/kraken2_db \
    --sylph_db /path/to/sylph.syldb \
    --diamond_db /path/to/amr.dmnd \
    --outdir ./results \
    -resume
