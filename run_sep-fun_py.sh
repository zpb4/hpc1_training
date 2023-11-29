
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=Pysep${SEED}\n\
#SBATCH -p normal\n\
#SBATCH --exclusive\n\
#SBATCH --export=ALL\n\
#SBATCH --nodes=1\n\
#SBATCH --output=Pysep${SEED}.txt\n\
#SBATCH --ntasks-per-node=80\n\

source ~/py-env/bin/activate\n\

python3 sep_experimentation.py"

echo -e $SLURM | sbatch
sleep 0.5
#done