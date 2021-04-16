#!/bin/bash -l

#SBATCH -A g2020012
#SBATCH --reservation g2020012_21
#SBATCH -p core -n 36
#SBATCH -t 5:00

module load gcc openmpi
# mpirun -np 4 bin/matmul /proj/g2020012/nobackup/matmul_indata/input3600.txt

for p in 1 2 4 9 16 25 36; do
    echo Num processes: $$p
    echo -n "$$p " >>$(SPEEDUP_OUT)
    for iter in $(seq 1 $(ITER)); do
        echo "Iteration: $$iter"
        mpirun --bind-to none -np $$p $(BINDIR)/$(BINFILE) $(SPEEDUP_IN) | xargs echo -ne "\t" >>$(SPEEDUP_OUT)
    done
    echo "" >>$(SPEEDUP_OUT)
done
