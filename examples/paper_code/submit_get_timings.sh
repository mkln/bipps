#!/bin/bash

#SBATCH --job-name simulation_get_exp
#SBATCH --output=logs/timings_%A_%a.out   # Standard output log (%A=job ID, %a=array index)
#SBATCH --error=logs/timings_%A_%a.err    # Standard error log (%A=job ID, %a=array index)
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00
#SBATCH --array=1-30
#SBATCH --mem=8G
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=END,FAIL

module load Rgeospatial/4.3.2-mkl-2024-01-15

Rscript examples/paper_code/get_timings.R $SLURM_ARRAY_TASK_ID
