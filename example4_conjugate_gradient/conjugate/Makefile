###############################################################################
# Conjugate gradient (CG) method
###############################################################################

# Project directories
BINDIR := bin
SRCDIR := src
TESTDIR := test

# C compliter
CC := mpicc


# Language standard
STD := -std=gnu99


# Stack protection
STACK := -fstack-protector-all -Wstack-protector


# Debug
DEBUG := -g3 -DDEBUG=1


# Dependency libraries 
LIBS = -lm

CFLAGS = -O3 $(STD) $(STACK)
BINFILE = conjugate


all: $(BINFILE)

$(BINFILE): $(SRCDIR)/$(BINFILE).c $(SRCDIR)/$(BINFILE).h
	$(CC) -o $(BINDIR)/$@ $< $(DEBUG) $(CFLAGS) $(LIBS)

test: $(BINFILE)  
	mpirun --oversubscribe -np 9 $(BINDIR)/$(BINFILE) 256

clean:
	$(RM) $(BINDIR)/*

