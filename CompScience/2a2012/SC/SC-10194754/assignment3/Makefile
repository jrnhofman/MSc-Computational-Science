CC=mpicc
CFLAGS=-O2 -Wall -o
LIB=-lapr

all: diff.o

diff.o: diff.c
	$(CC) $(CFLAGS) diff.o diff.c -L. $(LIB)

clean:
	rm -rf *.o all
