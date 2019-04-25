CC=mpicc
CFLAGS=-O2 -Wall -o
LIB=-lapr

all: mult.o

mult.o: mult.c
	$(CC) $(CFLAGS) mult.o mult.c

clean:
	rm -rf *.o all
