/* This is a template program for the third assignment: time-dependent diffusion. The space domain decompostion is using the strip-wise decomposition. You need to complete the missing codes(PART A,B,C,D,E,F). In many aspects you can learn from the vibrating string problem and maybe also the vector matrix product problem. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#define EXIT_SUCCESS 0
#endif
#define TRUE 1
#define FALSE 0
#define N 100           /*the number of grid data points in each row or column*/
#define dt 0.00001      /*the time stepsize*/
#define iter 1000       /* the number of iterations and it depends on both the simulation time and the time stepsize*/
static double *curval;  /* data grid point values for next time step*/
static double *newval;  /*current values of grid points on the concentration field*/
static MPI_Comm ring;   /*mpi communicator for ring/line topology */
static int nvalues;     /* number of values in array curval and newval  */
static int npoints;     /* number of grid data points for local process         */
static int nprocs;      /* number of processes/mpi-tasks in spmd implementation */
static int rank;        /* rank of this mpi process in ring/line topology       */
static int left;        /* rank of left mpi process in ring/line topology       */
static int right;        /* rank of right mpi process in ring/line topology      */
static double result[N]; /*this array is to store the elements of a column which you take from the computed concentration field*/

static void
compute(void)  /*
{
    
    PART   F

} 


int
main(int argc, char *argv[])
{
    int dims[1], periods[1], coords[1], neighbour[1];
    int i,j,k, offset;
    MPI_Status status;

    /* initialize mpi,get the environment and create ring topology */

 
    PART   A


    /* compute data decomposition and dynamic memory allocation of arrays*/


    PART   B


    /*Initialization of array elements */

    PART   C


    for(i=0;i<iter;i++)   compute();


    /*Take a column from the computed concentration field and then collect the result into processor 0 and print out*/

    PART   D

      /*Free the allocated memory and terminate MPI*/
    

    PART   E

    exit(EXIT_SUCCESS);
}






