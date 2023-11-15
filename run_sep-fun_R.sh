
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=Rsep${SEED}\n\
#SBATCH -p normal\n\
#SBATCH --export=ALL\n\
#SBATCH --nodes=1\n\
#SBATCH --exclusive\n\
#SBATCH --output=Rsep${SEED}.txt\n\
#SBATCH --ntasks-per-node=80\n\

module load R\n\

Rscript sep_experimentation.R"

echo -e $SLURM | sbatch 
sleep 0.5
