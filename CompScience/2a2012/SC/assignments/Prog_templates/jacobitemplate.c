/* This is the template program for the fourth assignment: the Jacobi iteration. You need to complete the missing codes(PART A,B,C,D,E,F,G).In many aspects you can learn from your program for the time dependent diffusion (i.e., the third assignment).*/

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
#define N 20             /*the number of grid data points in each row or column*/
#define sigma 1e-6       /*the stopping tolerance*/
#define niter 4
static double *curval;   /*the estimated concentration values of grid points in the current iteration */
static double *newval;   /*the estimated concentration values of grid points in the next iteration */
static MPI_Comm ring;    /*mpi communicator for ring/line topology */
static int nvalues;      /* number of rows in array curval and newval  */
static int npoints;      /* number of grid data rows for local process         */
static int nprocs;       /* number of processes/mpi-tasks in spmd implementation */
static int rank;         /* rank of this mpi process in ring/line topology       */
static int left;         /* rank of left mpi process in ring/line topology       */
static int right;        /* rank of right mpi process in ring/line topology      */
static double result[N];  /*this array is to store the elements of a column which you take from the computed concentration field*/
static double delta, maxdif;  /*delta is the local maximum difference between two iterations, maxdif is the global maximum difference*/

static void
compute(void)
{

    PART   G

} 


int
main(int argc, char *argv[])
{
    int dims[1], periods[1], coords[1], neighbour[1];
    int i,j,k, offset;
    long  iter;        /*the number of iterations*/
    MPI_Status status;


    /* initialize mpi,get the environment and create ring topology */

 
    PART   A


    /* compute data decomposition and dynamic memory allocation of arrays*/


    PART   B


    /*Initialization of array elements */

    PART   C
 

      iter=0;


    /*Implement the stopping criterion. You should use MPI_Allreduce to calculate the global maximum difference between two iterations*/
    

    PART   D

    /*Print out the number of iterations. Take a column from the computed concentration field and then collect the result into processor 0 and print out*/


    PART   E

      /*Free the allocated memory and terminate MPI*/
    

    PART   F

    exit(EXIT_SUCCESS);
}



