#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>

#include "apr.h"

struct parm parms;

int getparms(int argc, char *argv[])
{
	int rank = -1;
	int ret = -1;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (rank == 0)
	{
		ret = parseargs(argc, argv);
	}
	
	MPI_Bcast(&ret, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (ret != 0)
		return ret; // all procs propagate error from parseargs
	
	// quick 'n' dirty of broadcasting
	
	// send parms struct to all processes
	if (MPI_SUCCESS != MPI_Bcast((double*)(&parms.stime), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD))
		return 2;
	if (MPI_SUCCESS != MPI_Bcast((double*)(&parms.delta), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD))
		return 2;
	if (MPI_SUCCESS != MPI_Bcast((int*)(&parms.ntotal), 1, MPI_INT, 0, MPI_COMM_WORLD))
		return 2;
	if (MPI_SUCCESS != MPI_Bcast((int*)(&parms.freq), 1, MPI_INT, 0, MPI_COMM_WORLD))
		return 2;
	
	return 0;
}

/* Parse the commandline arguments */
int parseargs(int argc, char* argv[])
{
    char c;

    //Initialize valuese
    parms.ntotal = 0;
    parms.delta = 0.0;
    parms.stime = 0;
    
    //Read commandline
    while((c = getopt(argc, argv, "n:d:t:v:")) != -1) {
        switch (c) {
            case 'n':
                parms.ntotal = atoi(optarg);
                break;
            case 'd':
                parms.delta = atof(optarg);
                break;
            case 't':
                parms.stime = atof(optarg);
                break;
            case 'v':
                parms.freq = atoi(optarg);
                break;
        }
    }
   
    //Check if all required options are filled in
    if(!parms.delta || !parms.ntotal || !parms.stime) {
        printf("Usage: %s -n n -d d -t t -v v\n", argv[0]);
        return 1;       
    }
    return 0;
}



