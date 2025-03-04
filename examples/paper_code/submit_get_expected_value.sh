#!/bin/bash

#SBATCH --job-name simulation_get_exp
#SBATCH --output=logs/get_exp.out   # Standard output log (%A=job ID, %a=array index)
#SBATCH --error=logs/get_exp.err    # Standard error log (%A=job ID, %a=array index)
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00
#SBATCH --mem=8G
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=END,FAIL

module load Rgeospatial/4.3.2-mkl-2024-01-15

Rscript examples/get_expected_value.R
