#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH --reservation g2020012_21
#SBATCH -p core -n 36
#SBATCH -t 5:00

module load gcc openmpi
mpirun -np 4 bin/matmul /proj/g2020012/nobackup/matmul_indata/input3600.txt
