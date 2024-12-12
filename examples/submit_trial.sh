#!/bin/bash

#SBATCH --job-name simulation_trial
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=END,FAIL

module load Rgeospatial/4.3.2-mkl-2024-01-15

Rscript examples/simulation_runner.R
