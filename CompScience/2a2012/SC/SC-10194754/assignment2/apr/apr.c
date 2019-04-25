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
	
	// make sure conversion below is accurate (first bcast)
	assert(sizeof(enum { NONE=0, SINE=1, PLUCKED=2 }) == sizeof(int));
	
	// send parms struct to all processes
	if (MPI_SUCCESS != MPI_Bcast((int*)(&parms.method), 1, MPI_LONG, 0, MPI_COMM_WORLD))
		return 2;
	if (MPI_SUCCESS != MPI_Bcast((int*)(&parms.periods), 1, MPI_INT, 0, MPI_COMM_WORLD))
		return 2;
	if (MPI_SUCCESS != MPI_Bcast((double*)(&parms.location), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD))
		return 2;
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
    parms.freq = 0;
    parms.method = NONE;
    parms.ntotal = 0;
    parms.delta = 0.0;
    parms.stime = 0;
    
    //Read commandline
    while((c = getopt(argc, argv, "s:p:n:d:t:v:")) != -1) {
        switch (c) {
            case 's':
                parms.method = SINE;
                parms.periods = atoi(optarg);
                break;
            case 'p':
                parms.method = PLUCKED;
                parms.location = atof(optarg);
                break;
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
    //Get optimal tau if delta is not filled in
    if(!parms.delta) {
        parms.delta = 1.0 / (double)(parms.ntotal - 1);
        printf("Setting delta to %f\n", parms.delta);
    }
    
    //Check if all required options are filled in
    if(parms.method == NONE || !parms.ntotal || !parms.stime) {
        printf("Usage: %s {-s n | -p p} -n n -d d -t t [-v v]\n", argv[0]);
        return 1;       
    }
    return 0;
}

/* Initialization of the visualization procedure. 
 * This will open a file for the visualization 
 */
FILE *visualize_init()
{
	int rank = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    FILE *fp;
    char file[256];
    
    snprintf(file, 256, "string.m%d", rank);    
    fp = fopen(file, "w");
    
    return fp;
}

/*
 * Visualize the values of an iteration
 */
void visualize (FILE *fp, int  nvalues, double *values, float time)
{
    int rank = -1;
    int i;

	assert(nvalues >= 0);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Print Header
    if(rank == 0)
        fprintf(fp, "%f ", time);

    //Print Data
    for(i = 0; i < nvalues; i++)
        fprintf(fp, "%f ", values[i]);
    fprintf(fp, "\n");

    fflush(fp);
}

/* Combine all the files generated by the different processes */
void visualize_combine()
{
	MPI_Barrier(MPI_COMM_WORLD);

	int rank = -1;
	int nprocs = -1;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	//printf("debug: my rank in visualize_combine is %i, %i\n", rank, nprocs);

	if (rank == 0) // only one process needs to do it
	{
		int p;
		char command[4096];
		char buffer[4096];
    
		system("mv string.data string.data.bkp");

		//Concatenate data
		snprintf(command, 4096, "paste -d \"\"");

		//printf("debug: executing command 1:\n%s\n", command);

		for(p = 0; p < nprocs; p++)
		{
			strcpy(buffer, command);
			snprintf(command, 4096, "%s string.m%d", buffer, p);
		}

		//printf("debug: executing command 2:\n%s\n", command);

		strcpy(buffer, command);        
		snprintf(command, 4096, "%s > string.data", buffer);
		//snprintf(command, 4096, "%s", buffer);

		//printf("debug: executing command 3:\n%s\n", command);

		system(command);
            
		//Remove temporary files
		for(p = 0; p < nprocs; p++) {
			snprintf(command, 4096, "mv string.m%d string.m%d.bkp", p, p);
			system(command);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

/* Fit a line ax+b to the given data points and return the residual sum (sum of
 *  squared errors) */
double fit_line(int n, double *x, double *y, double *a, double *b)
{
	double SUMx, SUMy, SUMxy, SUMxx, SUMres, res, y_estimate ;
	int i;

	*a = 0.0;
	*b = 0.0;
	
	SUMx = 0; SUMy = 0; SUMxy = 0; SUMxx = 0;
	for (i=0; i<n; i++) {
		SUMx = SUMx + x[i];
		SUMy = SUMy + y[i];
		SUMxy = SUMxy + x[i]*y[i];
		SUMxx = SUMxx + x[i]*x[i];
	}
	*a = ( SUMx*SUMy - n*SUMxy ) / ( SUMx*SUMx - n*SUMxx );
	*b = ( SUMy - (*a)*SUMx ) / n;	
	
	SUMres = 0;
	for (i=0; i<n; i++) {
		y_estimate = (*a)*x[i] + *b;
		res = y[i] - y_estimate;
		SUMres = SUMres + res*res;
	}	
	
	return SUMres;
}