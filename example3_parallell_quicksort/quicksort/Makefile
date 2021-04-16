###############################################################################
# Parallel quicksort
###############################################################################

##############################
# Project directories
##############################
BINDIR := bin
SRCDIR := src
TESTDIR := test


##############################
# C compliter
##############################
CC := mpicc

##############################
# Language standard
##############################
STD := -std=gnu99

##############################
# Stack protection
##############################
STACK := -fstack-protector-all -Wstack-protector

##############################
# Debug
##############################
DEBUG := -g3 -DDEBUG=1

##############################
# Dependency libraries 
##############################
LIBS = -lm

CFLAGS = -O3 $(STD) $(STACK)
BINFILE = quicksort


TEST_SIZE = 4

#############################
# Speedup params
###########################

# Number of interations (for avering)
ITER := 10

# Sequence distrubution (0,1,2)
SEQ := 0

# 1D array size
LEN := 4000000

# Pivot method (0,1,2)
PIV := 0

# output file
SPEEDUP_OUT := test/speedup.txt

#############################
# Strong scale test
###########################
STRONG_LEN := 400000000
STRONG_OUT := test/strongscale.txt 

#############################
# Weak scale test
###########################
WEAK_OUT := test/weakscale.txt 

all: $(BINFILE)

$(BINFILE): $(SRCDIR)/$(BINFILE).c $(SRCDIR)/$(BINFILE).h
	$(CC) -o $(BINDIR)/$@ $< $(DEBUG) $(CFLAGS) $(LIBS)

test: $(BINFILE)  
	mpirun -np 8 bin/quicksort 2 400 1


speedup: $(BINFILE)
	for p in 1 4 16 32; do \
		echo Num processes: $$p ; \
		echo -n "$$p " >> $(SPEEDUP_OUT)  ; \
		for iter in `seq 1 $(ITER)`; do \
			echo iteration: $$iter; \
			mpirun --bind-to none -np $$p $(BINDIR)/$(BINFILE)  $(SEQ) $(LEN) $(PIV)| xargs echo -ne "\t">> $(SPEEDUP_OUT)  ; \
		done; \
		echo "" >> $(SPEEDUP_OUT); \
	done

strongscale: $(BINFILE)
	for p in 1 4 16 32; do \
		echo Num processes: $$p ; \
		echo -n "$$p " >> $(STRONG_OUT)  ; \
		for iter in `seq 1 $(ITER)`; do \
			echo iteration: $$iter; \
			mpirun --bind-to none -np $$p $(BINDIR)/$(BINFILE)  $(SEQ) $(LEN) $(PIV)| xargs echo -ne "\t">> $(STRONG_OUT)  ; \
		done; \
		echo "" >> $(STRONG_OUT); \
	done

weakscale: $(BINFILE)
	for p in 1 4 16 32; do \
		len=$$(( $$p * 10000000 )); \
		echo Num processes: $$p ; \
		echo -n "$$p " >> $(WEAK_OUT)  ; \
		for iter in `seq 1 $(ITER)`; do \
			echo iteration: $$iter; \
			mpirun --bind-to none -np $$p $(BINDIR)/$(BINFILE)  $(SEQ) $$len $(PIV)| xargs echo -ne "\t">> $(WEAK_OUT)  ; \
		done; \
		echo "" >> $(WEAK_OUT); \
	done

clean:
	$(RM) $(BINDIR)/*
