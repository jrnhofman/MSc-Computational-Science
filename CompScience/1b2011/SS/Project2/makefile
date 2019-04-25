CC = gcc
LIBS = -lgsl -lgslcblas -lm
OPT = -O3
CCFLAGS = $(LIBS)

all: main 

main: rng_test.c MT/mt19937ar.c
	$(CC) $(OPT) -o rng_test rng_test.c MT/mt19937ar.c $(CCFLAGS)

clean: main 

