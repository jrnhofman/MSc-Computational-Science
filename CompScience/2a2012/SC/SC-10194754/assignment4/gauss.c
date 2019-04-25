#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <sys/time.h>
#include "apr.h"

#define PI 3.14159265358979323846

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
    double old;
    localdiff = 0;
    int nhalve = (N-2)/2;
    int offset = nhalve*((rank+iter)%2);
    int offset2 = nhalve*((rank+1+iter)%2);

    /* Pass global array, which are copies of the array in main() */
    double (*curval)[nvalues][N] = (double (*)[nvalues][N]) curvalx;

    /* Apply finite difference to 'black' tiles and copy columns */
    if(iter%10==0)
    {
	for(i=1; i<npoints+1; i++)
	{
	    for(j=1+offset; j<nhalve+1+offset; j++)
	    {
		old = (*curval)[i][j];
		(*curval)[i][j] = 0.25*((*curval)[i+1][j] + (*curval)[i-1][j] + (*curval)[i][j-1] + (*curval)[i][j+1]);
		localdiff = fmax(localdiff,(*curval)[i][j]-old);
	    }
	}
    }
    else
    {
	for(i=1; i<npoints+1; i++)
	{
	    for(j=1+offset; j<nhalve+1+offset; j++)
	    {
		(*curval)[i][j] = 0.25*((*curval)[i+1][j] + (*curval)[i-1][j] + (*curval)[i][j-1] + (*curval)[i][j+1]);
	    }
	}
    }

    if(offset==0)
    {
	for(i=1; i<npoints+1; i++)
	{
	    (*curval)[i][N-1] = (*curval)[i][1];
	}
    }
    else
    {
	for(i=1; i<npoints+1; i++)
	{
	    (*curval)[i][0] = (*curval)[i][N-2];
	}
    }
	
    /* Send and receive boundary points */
    if(rank!=0)
    {
	MPI_Isend(&(*curval)[1][offset],nhalve+1,MPI_DOUBLE,up,0,line,&req[0]);
	MPI_Irecv(&(*curval)[0][offset2],nhalve+1,MPI_DOUBLE,up,0,line,&req[2]);
    }

    if(rank!=nprocs-1)
    {
	MPI_Isend(&(*curval)[npoints][offset],nhalve+1,MPI_DOUBLE,down,0,line,&req[1]);
	MPI_Irecv(&(*curval)[npoints+1][offset2],nhalve+1,MPI_DOUBLE,down,0,line,&req[3]);
    }

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

    /* for(i=0; i<npoints+2; i++) */
    /* { */
    /* 	for(j=0; j<N; j++) */
    /* 	{ */
    /* 	    printf("(%d) %f",rank,(*curval)[i][j]); */
    /* 	} */
    /* 	printf("\n"); */
    /* } */
    /* printf("\n"); */


    /* Apply finite difference to 'red' tiles and copy columns */
    if(iter%10==0)
    {
	for(i=1; i<npoints+1; i++)
	{
	    for(j=1+offset2; j<nhalve+1+offset2; j++)
	    {
		old = (*curval)[i][j];
		(*curval)[i][j] = 0.25*((*curval)[i+1][j] + (*curval)[i-1][j] + (*curval)[i][j-1] + (*curval)[i][j+1]);
		localdiff = fmax(localdiff,(*curval)[i][j]-old);
	    }
	}
    }
    else
    {
	for(i=1; i<npoints+1; i++)
	{
	    for(j=1+offset2; j<nhalve+1+offset2; j++)
	    {
		(*curval)[i][j] = 0.25*((*curval)[i+1][j] + (*curval)[i-1][j] + (*curval)[i][j-1] + (*curval)[i][j+1]);
	    }
	}
    } 

    if(offset2==0)
    {
	for(i=1; i<npoints+1; i++)
	{
	    (*curval)[i][N-1] = (*curval)[i][1];
	}
    }
    else
    {
	for(i=1; i<npoints+1; i++)
	{
	    (*curval)[i][0] = (*curval)[i][N-2];
	}
    }
	
    /* Send and receive boundary points */
    if(rank!=0)
    {
	MPI_Isend(&(*curval)[1][offset2],nhalve+1,MPI_DOUBLE,up,0,line,&req[0]);
	MPI_Irecv(&(*curval)[0][offset],nhalve+1,MPI_DOUBLE,up,0,line,&req[2]);
    }

    if(rank!=nprocs-1)
    {
	MPI_Isend(&(*curval)[npoints][offset2],nhalve+1,MPI_DOUBLE,down,0,line,&req[1]);
	MPI_Irecv(&(*curval)[npoints+1][offset],nhalve+1,MPI_DOUBLE,down,0,line,&req[3]);
    }

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


    /* for(i=0; i<npoints+2; i++) */
    /* { */
    /* 	for(j=0; j<N; j++) */
    /* 	{ */
    /* 	    printf("(%d) %f",rank,(*curval)[i][j]); */
    /* 	} */
    /* 	printf("\n"); */
    /* } */
    /* printf("\n"); */

} 


void visualize(int iter)
{

    /* initialize file for data */
    FILE *fp;
    char file[256];
    snprintf(file, 256, "data.%dgauss%.0e",parms.ntotal,parms.threshold);    
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

    /* linking the arrays above with the global ones */
    curvalx = curval;

    /* Initialization of grid, set first row to 1, rest to 0 */
    for(i=0; i<nvalues; i++)
    {
	for(j=0; j<N; j++)
	{
	    (*curval)[i][j] = 0;
	}
    }

    if(rank==0)
    {
	for(j=0; j<N; j++)
	{
	    (*curval)[0][j] = 1;
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
    //int niters = 1;
    //for(iter=0; iter<niters; iter++) {
    while(globaldiff>parms.threshold)
    {
         /* visualize string at the first iteration and then every t.freq
         * time steps.
         */
        compute();
	
        /* Reduce localdiff to globaldiff */
	if(iter%10==0)
	{
	    MPI_Allreduce(&localdiff,&globaldiff,1,MPI_DOUBLE,MPI_MAX,line);
	}

	iter++;
    }

    /* Final visualization */
    if(parms.freq==1)
	visualize(iter);
 	
    /* Record elapsed time */
    if(rank==0)
    {
	printf("Iteration: %d, Time: %f, Threshold: %f\n",iter-1,(iter-1)*parms.delta,parms.threshold);
	gettimeofday(&endTime,NULL);
	double time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	printf("Matrix size: %d\t Number of processes: %d\t Time per iteration: %f s\n",parms.ntotal,nprocs,time/(iter-1));
    }

    /* Cleaning up */
    free(curval);
    free(req);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}


