CC = g++
LIBS = -lm
OPT = -O3
CCFLAGS = $(LIBS)

all: gen sim

sim: TSP_sim.cpp SFMT/SFMT.c
	$(CC) $(OPT) -o TSP_sim TSP_sim.cpp SFMT/SFMT.c $(CCFLAGS)

gen: TSP_gen.cpp SFMT/SFMT.c
	$(CC) $(OPT) -o TSP_gen TSP_gen.cpp SFMT/SFMT.c $(CCFLAGS)

clean: gen sim

