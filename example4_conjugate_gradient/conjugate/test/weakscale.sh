#!/bin/bash -l


module load gcc openmpi

n=1024
BINFILE=../bin/conjugate

for p in 1 2 3 4 5 6 
do
	echo num procs: $(( $p*$p )) >> weakscale.txt
	for iter in {1..10}  
	do
		mpirun --oversubscribe -np $(( $p*$p )) $BINFILE $(( $n*p )) >> weakscale.txt
	done
done
