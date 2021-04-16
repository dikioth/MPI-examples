###############################################################################
# Makefile for assignment 2 
# Parallel and Distributed Computing 2020.
###############################################################################

# C compliter
CC := mpicc

# Language standard
STD := -std=c99

# Stack protection
STACK := -fstack-protector-all -Wstack-protector

# Debug
DEBUG := -g3 -DDEBUG=1

# Dependency libraries 
LIBS = -lm

CFLAGS = -O3 $(STD) $(STACK)

BINFILE = matmul

# Project directories
BINDIR := bin
SRCDIR := src
TESTDIR := test

TEST_SIZE = 4

# Speedup params
ITER := 10
SPEEDUP_IN := ../test_data/matmul_indata/input3600.txt
SPEEDUP_OUT := test/speedup.txt


all: $(BINFILE)

$(BINFILE): $(SRCDIR)/matmul.c $(SRCDIR)/matmul.h
	$(CC) -o $(BINDIR)/$@ $< $(DEBUG) $(CFLAGS) $(LIBS)

test: $(BINFILE)  
	$(BINDIR)/$(BINFILE) $(TESTDIR)/input4.txt $(TESTDIR)/temp1.txt
	diff  -w  $(TESTDIR)/temp1.txt $(TESTDIR)/output4_ref.txt
	mpirun -np 4 $(BINDIR)/$(BINFILE) $(TESTDIR)/input4.txt $(TESTDIR)/temp2.txt
	diff  -w  $(TESTDIR)/temp2.txt $(TESTDIR)/output4_ref.txt
	rm $(TESTDIR)/temp*.txt 
	echo "test passed!"

compare: 
	mpirun -np $$p $(BINDIR)/$(BINFILE) $$ifile $(TESTDIR)/temp.txt
	diff  -w  $(TESTDIR)/temp.txt $$ofile
	rm $(TESTDIR)/temp.txt

speedup: $(BINFILE)
	for p in 1 4 9 16; do \
		echo Num processes: $$p ; \
		echo -n "$$p " >> $(SPEEDUP_OUT)  ; \
		for iter in `seq 1 $(ITER)`; do \
			echo iteration: $$iter; \
			mpirun --bind-to none -np $$p $(BINDIR)/$(BINFILE) $(SPEEDUP_IN) | xargs echo -ne "\t">> $(SPEEDUP_OUT)  ; \
		done; \
		echo "" >> $(SPEEDUP_OUT); \
	done

clean:
	$(RM) $(BINDIR)/*
