#!/bin/bash

#SBATCH -t 100:00:00
#SBATCH --job-name=par-ex
#SBATCH -p normal
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --output=par-ex.txt
#SBATCH --ntasks-per-node=80

module load R/4.1.2

Rscript parallel_examples.R


