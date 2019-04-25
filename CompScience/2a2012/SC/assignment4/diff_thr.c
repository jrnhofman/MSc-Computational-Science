#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <sys/time.h>
#include "apr.h"

#define PI 3.14159265358979323846

void *newvalx;
void *curvalx;
static MPI_Comm line;  /* mpi communicator for line topology       */
static int nvalues; /* number of values in array curval, oldval and newval  */
static int npoints; /* number of grid data points for local process         */
static int nprocs;  /* number of processes/mpi-tasks in spmd implementation */
static int rank;    /* rank of this mpi process in ring/line topology       */
static int up;    /* rank of upper mpi process in ring/line topology       */
static int down;   /* rank of down mpi process in ring/line topology      */
static int N; /* matrix size */

double globaldiff;
double localdiff;
int iter;

static MPI_Request *req;
double *swap;

struct parm parms; /* struct with values from command line */ 

FILE   *fp;         /* file pointer for visualization                       */

void compute(void)
{
    int i,j;
    localdiff = 0;
    double D = parms.delta * parms.ntotal * parms.ntotal; /*D dt / dx^2*/

    /* Pass global arrays, which are copies of the arrays in main() */
    double (*newval)[nvalues][N] = (double (*)[nvalues][N]) newvalx;
    double (*curval)[nvalues][N] = (double (*)[nvalues][N]) curvalx;

    /* Apply finite difference to border rows and copy to halo columns */
    for(j=1; j<N-1; j++)
    {
	(*newval)[1][j] = (*curval)[1][j] + D * 
	    ((*curval)[2][j] + (*curval)[0][j] + (*curval)[1][j-1] + (*curval)[1][j+1] - 4*(*curval)[1][j]);
	(*newval)[npoints][j] = (*curval)[npoints][j] + D * 
	    ((*curval)[npoints+1][j] + (*curval)[npoints-1][j] + (*curval)[npoints][j-1] + (*curval)[npoints][j+1] - 4*(*curval)[npoints][j]);
    }
    (*newval)[1][0] = (*newval)[1][N-2];
    (*newval)[1][N-1] = (*newval)[1][1];
    (*newval)[npoints][0] = (*newval)[npoints][N-2];
    (*newval)[npoints][N-1] = (*newval)[npoints][1];

    /* Send and receive boundary points */
    if(rank!=0)
    {
	MPI_Isend(&(*newval)[1][0],N,MPI_DOUBLE,up,0,line,&req[0]);
	MPI_Irecv(&(*newval)[0][0],N,MPI_DOUBLE,up,0,line,&req[2]);
    }

    if(rank!=nprocs-1)
    {
	MPI_Isend(&(*newval)[npoints][0],N,MPI_DOUBLE,down,0,line,&req[1]);
	MPI_Irecv(&(*newval)[npoints+1][0],N,MPI_DOUBLE,down,0,line,&req[3]);
    }
    
    /* Apply finite difference to rest of rows and copy to halo columns */
    for(i=2; i<npoints; i++)
    {
	for(j=1; j<N-1; j++)
	{
	    (*newval)[i][j] = (*curval)[i][j] + D * 
		((*curval)[i+1][j] + (*curval)[i-1][j] + (*curval)[i][j-1] + (*curval)[i][j+1] - 4*(*curval)[i][j]);
	}
    }

    for(i=2; i<npoints; i++)
    {
	(*newval)[i][0] = (*newval)[i][N-2];
	(*newval)[i][N-1] = (*newval)[i][1];
    }

    /* Calculate local maximum difference */
    if(iter%10==0)
    {
	for(i=1; i<npoints+1; i++)
	{
	    for(j=1; j<N-1; j++)
	    {
		localdiff = fmax(localdiff,(*newval)[i][j]-(*curval)[i][j]);
	    }
	}
    }

    /* for(i=1;i<npoints+1;i++) */
    /* { */
    /* 	for(j=1;j<N-1;j++) */
    /* 	{ */
    /* 	    printf("%f ",(*curval)[i][j]); */
    /* 	} */
    /* 	printf("\n"); */
    /* } */
    /* printf("\n"); */
    /* for(i=1;i<npoints+1;i++) */
    /* { */
    /* 	for(j=1;j<N-1;j++) */
    /* 	{ */
    /* 	    printf("%f ",(*newval)[i][j]); */
    /* 	} */
    /* 	printf("\n"); */
    /* } */
    /* printf("\n"); */

    /* Wait for communication to be completed */
    if(rank!=0)
    {
	MPI_Wait(&req[0],MPI_STATUS_IGNORE);
	MPI_Wait(&req[2],MPI_STATUS_IGNORE);
    }
    if(rank!=nprocs-1)
    {
	MPI_Wait(&req[1],MPI_STATUS_IGNORE);
	MPI_Wait(&req[3],MPI_STATUS_IGNORE);
    }

    /* Swap arrays */
    swap = curvalx;
    curvalx = newvalx;
    newvalx = swap;
} 


void visualize(int iter)
{

    /* initialize file for data */
    FILE *fp;
    char file[256];
    snprintf(file, 256, "data.i%dp%d",iter,rank);    
    fp = fopen(file, "w");

    double (*curval)[nvalues][N] = (double (*)[nvalues][N]) curvalx;

    /* Print current iteration in file */
    int i, j;
    for(i=1; i<nvalues-1; i++)
    {
	for(j=1; j<N-1; j++)
	{
	    fprintf(fp,"%f ",(*curval)[i][j]);
	}
	fprintf(fp,"\n");
    }

    fclose(fp);
}

int main(int argc, char *argv[])
{
    int dims[1], periods[1];
    int i, j;

    req = (MPI_Request*) malloc(sizeof(MPI_Request)*4);

    /* initialize mpi */
    MPI_Init(&argc, &argv);

    /* rank process 0 in world communicator receives parameters on
     * command line, parses them and pass them to other processes.
     * the other mpi processes receive the parameters.*/
    getparms(argc,argv);
   
    /* get number of mpi processes (nprocs) and initialize ring/line
     * topology (the communicator stored in ring variable).  Find
     * out own rank and left and right processes ranks on ring.
     */
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    dims[0] = nprocs; 
    periods[0] = 0;
    MPI_Cart_create(MPI_COMM_WORLD,1,dims,periods,1,&line);
    MPI_Comm_rank(line,&rank);
    MPI_Cart_shift(line,1,1,&up,&down);

    MPI_Barrier(line);

    /* Start timer */
    struct timeval startTime;
    struct timeval endTime;
    if(rank==0)
    {
    	gettimeofday(&startTime,NULL);
    }

    /* compute data decomposition;
     * - npoints should hold the number of grid data points on the string
     *   for the local computing node.
     * - nvalues is the number of elements in the array (which is npoints
     *   plus the neighbourhood.
     */

    N = parms.ntotal+2; // total number of points
    npoints = parms.ntotal/nprocs+(parms.ntotal%nprocs+nprocs-rank-1)/nprocs;
    nvalues = npoints + 2;

    /* allocation of dataset */
    double (* curval)[nvalues][N] =				\
	(double (*)[nvalues][N]) malloc(nvalues*(N)*sizeof(double));
    double (* newval)[nvalues][N] =				\
	(double (*)[nvalues][N]) malloc(nvalues*(N)*sizeof(double));

    /* linking the arrays above with the global ones */
    curvalx = curval;
    newvalx = newval;

    /* Initialization of grid, set first row to 1, rest to 0 */
    for(i=0; i<nvalues; i++)
    {
	for(j=0; j<N; j++)
	{
	    (*curval)[i][j] = 0;
	    (*newval)[i][j] = 0;
	}
    }

    if(rank==0)
    {
	for(j=0; j<N; j++)
	{
	    (*curval)[0][j] = 1;
	    (*newval)[0][j] = 1;
	}
    }

    /* computation loop, determine number of iteration based on simulation
     * time and time advance per iteration.  call visualization at determined
     * number of iterations and when loop is done.  don't call visualization
     * if not requested (parameter is equal to zero).
     */

    globaldiff = parms.threshold+1;
    iter = 0;
    //niters = (int) ceil((double) parms.stime / parms.delta);
    //int niters = 1000;
    //for(iter=0; iter<niters; iter++) {
    while(globaldiff>parms.threshold)
    {
         /* visualize string at the first iteration and then every t.freq
         * time steps.
         */
        if(parms.freq > 0 && iter % parms.freq == 0)
	    visualize(iter); 
        compute();
	
        /* Reduce localdiff to globaldiff */
	if(iter%10==0)
	{
	    MPI_Allreduce(&localdiff,&globaldiff,1,MPI_DOUBLE,MPI_MAX,line);
	}

	iter++;
    }

    /* Final visualization */
    if(parms.freq > 0)
	visualize(iter);
 	
    /* Record elapsed time */
    if(rank==0)
    {
	//printf("Iteration: %d, Time: %f, Threshold: %f\n",iter-1,(iter-1)*parms.delta,parms.threshold);
	gettimeofday(&endTime,NULL);
	double time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	printf("Matrix size: %d\t Number of processes: %d\t Time per iteration: %f s\n",parms.ntotal,nprocs,time/(iter-1));
    }

    /* Cleaning up */
    free(curval);
    free(newval);
    free(req);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}


