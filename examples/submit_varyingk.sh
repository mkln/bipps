#!/bin/bash

#SBATCH --job-name simulation_varyingk
#SBATCH --array=1-30
#SBATCH --output=logs/job_%A_%a.out   # Standard output log (%A=job ID, %a=array index)
#SBATCH --error=logs/job_%A_%a.err    # Standard error log (%A=job ID, %a=array index)
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=END,FAIL

module load Rgeospatial/4.3.2-mkl-2024-01-15

Rscript examples/simulations_varyingk_slurm_runner.R ${SLURM_ARRAY_TASK_ID}
