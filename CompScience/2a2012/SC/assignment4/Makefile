CC=mpicc
CFLAGS=-O2 -Wall -o
LIB=-lapr

all: diff_thr.o jacobi.o gauss.o

diff_thr.o: diff_thr.c
	$(CC) $(CFLAGS) diff_thr.o diff_thr.c -L. $(LIB)
jacobi.o: jacobi.c
	$(CC) $(CFLAGS) jacobi.o jacobi.c -L. $(LIB)
gauss.o: gauss.c
	$(CC) $(CFLAGS) gauss.o gauss.c -L. $(LIB)

clean:
	rm -rf *.o all
