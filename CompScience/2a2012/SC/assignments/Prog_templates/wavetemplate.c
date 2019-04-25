#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <apr.h>

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

struct parm parms; /* structure with parameters which describe initial
                    * string configuration and computation parameters;
                    * see further apr.h file for fields.
                    */

FILE   *fp;         /* file pointer for visualization                       */

static void
compute(void)
{
    
    double tau2; /* squared tau value, where tau equals c * dt / dx */
    tau2 = 1.0 * parms.delta / 1.0;
    tau2 = tau2 * tau2;
   
   PART D
} 


int
main(int argc, char *argv[])
{
    int dims[1], periods[1], coords[1], neighbour[1];
    int i, iter, niters, offset;

    /* initialize mpi */
    MPI_Init(&argc, &argv);
    /* rank process 0 in world communicator receives parameters on
     * command line, parses them and pass them to other processes.
     * the other mpi processes receive the parameters.*/
     
   PART A1
   
    /* get number of mpi processes (nprocs) and initialize ring/line
     * topology (the communicator stored in ring variable).  Find
     * out own rank and left and right processes ranks on ring.
     */
    
    

   PART A2

  

    /* compute data decomposition;
     * - npoints should hold the number of grid data points on the string
     *   for the local computing node.
     * - nvalues is the number of elements in the array (which is npoints
     *   plus the neighbourhood.
     * - offset is the starting position of the local grid data points
     *   in the global grid data points.
     */
    
    npoints = ...
    nvalues = ...
    offset=...

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
    
    PART C

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
	/* visualize end */
        visualize(fp, npoints, &curval[1], iter*parms.delta);

    /* finish visualization, deallocate memory and quit from mpi */
    if(parms.freq > 0)
        fclose(fp);

    free(curval);
    free(oldval);
    free(newval);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}





