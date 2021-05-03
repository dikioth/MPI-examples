#!/bin/bash -l

module load gcc openmpi

BINFILE=../bin/conjugate
n=2048
for p in 1 4 9 16 25 36 
do
	echo num procs: $p >> speedup.txt
	for iter in {1..10}  
	do
		mpirun --oversubscribe -np $p $BINFILE $n >> speedup.txt
	done
done
