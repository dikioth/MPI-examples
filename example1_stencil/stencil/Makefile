###############################################################################
# Stencil 
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm
BINDIR = bin
SRCDIR = src
BIN = stencil

# NUMBER ITERATIONS (FOR TIME AVERAGE)
ITER = 10

all: $(BIN)

$(BIN): $(SRCDIR)/stencil.c $(SRCDIR)/stencil.h
	$(CC) $(CFLAGS) -o $(BINDIR)/$@ $< $(LIBS)

test: $(BIN)  
	mpirun -np 4 $(BINDIR)/$(BIN)  test/input96.txt test/test1.txt 1 > /dev/null
	mpirun -np 4 $(BINDIR)/$(BIN)  test/input96.txt test/test2.txt 4 > /dev/null
	mpirun -np 4 $(BINDIR)/$(BIN)  test/input120.txt test/test3.txt 2 > /dev/null
	diff --ignore-all-space test/test1.txt test/output96_1_ref.txt
	diff --ignore-all-space test/test2.txt test/output96_4_ref.txt
	diff --ignore-all-space test/test3.txt test/output120_2_ref.txt
	rm test/test*.txt

clean:
	$(RM) $(BIN)
