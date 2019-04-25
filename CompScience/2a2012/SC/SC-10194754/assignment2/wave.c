#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <sys/time.h>
#include "apr.h"

#define PI 3.14159265358979323846

static double *curval; /* current values of grid points on string       */
static double *oldval; /* previous time step values of data grid points */
static double *newval; /* data grid point values for next time step     */
static MPI_Comm ring;  /* mpi communicator for ring/line topology       */
static int nvalues; /* number of values in array curval, oldval and newval  */
static int npoints; /* number of grid data points for local process         */
static int nprocs;  /* number of processes/mpi-tasks in spmd implementation */
static int rank;    /* rank of this mpi process in ring/line topology       */
static int left;    /* rank of left mpi process in ring/line topology       */
static int right;   /* rank of right mpi process in ring/line topology      */

static MPI_Request *req;
double *swap;

struct parm parms; /* structure with parameters which describe initial
                    * string configuration and computation parameters;
                    * see further apr.h file for fields.
                    */

FILE   *fp;         /* file pointer for visualization                       */

static void compute(void)
{
    int i;
    double tau2 = parms.delta * parms.ntotal; /*c^2 dt^2 / dx^2*/
    tau2 = tau2 * tau2;
   
    /* Apply finite difference */
    for(i=1; i<nvalues-1; i++)
    {
	newval[i] = 2*curval[i] - oldval[i] +
	    tau2 * (curval[i+1] - 2*curval[i] + curval[i-1]);
    }

    /* Move all arrays */
    swap = oldval;
    oldval = curval;
    curval = newval;
    newval = swap;

    /* Send and receive boundary points */
    MPI_Isend(&curval[1],1,MPI_DOUBLE,left,1,ring,&req[0]);
    MPI_Isend(&curval[npoints],1,MPI_DOUBLE,right,0,ring,&req[1]);
    MPI_Irecv(&curval[0],1,MPI_DOUBLE,left,0,ring,&req[2]);
    MPI_Irecv(&curval[npoints+1],1,MPI_DOUBLE,right,1,ring,&req[3]);

    /* Wait for communication to be completed */
    MPI_Waitall(4,req,MPI_STATUSES_IGNORE);
} 


int main(int argc, char *argv[])
{
    int dims[1], periods[1];
    int i, iter, niters, offset;
    double fraction, slope_up, slope_down, x;

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
    periods[0] = 1;
    MPI_Cart_create(MPI_COMM_WORLD,1,dims,periods,1,&ring);
    MPI_Comm_rank(ring,&rank);
    MPI_Cart_shift(ring,1,1,&left,&right);

    MPI_Barrier(ring);

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
     * - offset is the starting position of the local grid data points
     *   in the global grid data points.
     */
    
    npoints = parms.ntotal/nprocs+(parms.ntotal%nprocs+nprocs-rank-1)/nprocs;
    nvalues = npoints + 2;
    if(npoints>parms.ntotal/nprocs) { offset = npoints * rank; }
    else { offset = parms.ntotal - (nprocs - rank)*npoints; }

    /* allocation of dataset and initialize visualisation */
    if(!(curval = malloc(sizeof(double) * nvalues)) ||
       !(oldval = malloc(sizeof(double) * nvalues)) ||
       !(newval = malloc(sizeof(double) * nvalues))) {
        fprintf(stderr,"%s: out of memory\n",argv[0]);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    if(parms.freq > 0)
        /* if frequency is zero, don't create visualization file */
	fp = visualize_init(rank);

    /* initialization of string, initialize all grid data points in the
     * array according to the parameters set in parms. I.e.; either a
     * sine or a plucked string with eithers a number of periods or
     * plucked position.
     * use variables offset and parms.ntotal to determin the local x
     * position.
     */

    /* Sinus */
    if(parms.method==1)
    {
	fraction = (double)parms.periods*2*PI/parms.ntotal;
	for(i=0; i<nvalues; i++)
	{
	    curval[i] = sin((double)(offset+i-1)*fraction);
	}
    }

    /* Plucked */
    if(parms.method==2)
    {
	slope_up = 1./(parms.location);
	slope_down = -1./(1-parms.location);
	for(i=0; i<nvalues; i++)
	{
	    x = (double)(offset+i-1)/parms.ntotal;
	    if(x<=parms.location)
		curval[i] = x * slope_up;
	    else
		curval[i] = 1 + ((x-parms.location) * slope_down);
	}
	if(rank==0)
	{
	    x = (double)(parms.ntotal-1)/parms.ntotal;
	    curval[0] = 1 + ((x-parms.location) * slope_down);
	}
	else if(rank==nprocs-1)
	{
	    curval[nvalues-1] = 0;
	}
    }

    /* Set old values to current values for initialisation */
    memcpy(oldval,curval,sizeof(double)*nvalues);

    /* computation loop, determin number of iteration based on simulation
     * time and time advance per iteration.  call visualization at determined
     * number of iterations and when loop is done.  don't call visualization
     * if not requested (parameter is equal to zero).
     */

    niters = (int) ceil((double) parms.stime / parms.delta);
    for(iter=0; iter<niters; iter++) {
        /* visualize string at the first iteration and then every t.freq
         * time steps.
         */
        if(parms.freq > 0 && iter % parms.freq == 0)
            visualize(fp, npoints, &curval[1], iter*parms.delta);
        compute();
    }
    if(parms.freq > 0)
    {
	/* visualize end and combine */
        visualize(fp, npoints, &curval[1], iter*parms.delta);
	visualize_combine();
	fclose(fp);
    }
	
    /* Record elapsed time */
    if(rank==0)
    {
	gettimeofday(&endTime,NULL);
	double time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	printf("String size: %d\t Number of processes: %d\t Time: %f s\n",parms.ntotal,nprocs,time);
    }

    /* Cleaning up */
    free(curval);
    free(oldval);
    free(newval);
    free(req);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
