CC=mpicc
CFLAGS=-O2 -Wall -o
LIB=-lapr

all: wave.o

wave.o: wave.c
	$(CC) $(CFLAGS) wave.o wave.c -L. $(LIB)

clean:
	rm -rf *.o all
