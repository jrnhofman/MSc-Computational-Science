CC = g++
WARN = -Wall

CCFLAGS = $(WARN)

all: main graphics

main: main.cpp MT/mt19937ar.c
	$(CC) -o main main.cpp MT/mt19937ar.c $(CCFLAGS)

graphics: graphics.cpp
	$(CC) -o graphics graphics.cpp $(CCFLAGS)
	
clean: main graphics

