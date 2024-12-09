#!/bin/bash

#SBATCH --job-name simulation_trial
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=END,FAIL

module load Rgeospatial

Rscript examples/simulation_runner.R
