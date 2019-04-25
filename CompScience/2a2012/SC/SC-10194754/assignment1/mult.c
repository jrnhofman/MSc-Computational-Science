#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

int main(int argc, char *argv[])
{
    int i, j, rank;
    int p, n;  /* p is the number of processors */
    int nrows, offset;  /* nrows is the number of rows which the 
			     current node is responsible to calculate */

    if(argc==1) 
    {
	printf("Setting the matrix size to 100!\n");
	n = 100;
    }
    else { n = atoi(argv[1]); }

    struct timeval startTime;
    struct timeval endTime;

    /* Initialize MPI */
    MPI_Init(NULL,NULL);

    /* Get environment, i.e. the number of processors and the i.d. */
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Request *req = (MPI_Request*) malloc(sizeof(MPI_Request)*(p-1));
    if(p > n)
    {
	if(rank==0)
	    printf("Too many processes, maximum allowed: %d/n", n);
        MPI_Abort(MPI_COMM_WORLD,0);
    }

    /* Begin time measurement */
    if(rank==0)
    {
    	gettimeofday(&startTime,NULL);
    }

    /* Compute data decomposition, i.e., the values of nrows and offset*/
    nrows = n/p+(n%p+p-rank-1)/p;
    if(nrows>n/p) { offset=nrows*rank; }
    else { offset=n-(p-rank)*nrows; }
 
    /* Initialization of the elements of the vectors and the matrix */
    double *vector = (double*) malloc(sizeof(double)*n);
    double *res = (double*) malloc(sizeof(double)*n);
    double (*matrix)[nrows][n] = \
	(double (*)[nrows][n]) malloc(nrows*n*sizeof(double));
    for(i=offset;i<offset+nrows;i++) { vector[i] = i; }
    for(i=0;i<nrows;i++)
    {
	for(j=0;j<n;j++) 
	{
	    (*matrix)[i][j]=i+offset;
	}
    }
    
    /* Calculate length and displacement matrix */
    int *length = (int*) malloc(sizeof(int)*p);
    int *displ = (int*) malloc(sizeof(int)*p);
    for(i=0; i<p; i++)
    {
	length[i] = (n/p + (n%p+p-i-1)/p);
	if(i!=0) { displ[i] = displ[i-1] + length[i-1]; }
	else { displ[i] = 0; }
    }

    /* Gather all the data to make a complete vector on all processes*/
    MPI_Allgatherv(&vector[offset],nrows,MPI_DOUBLE,
		  vector,length,displ,MPI_DOUBLE,MPI_COMM_WORLD);

    /* /*Print the matrix and vector */ 
    /* if(rank==0) */
    /* { */
    /* 	printf("the elements of the vector are:\n"); */
    /* 	for(i=0;i<n;i++) */
    /* 	{  */
    /* 	    printf("%0.0f  ",vector[i]); */
    /* 	} */
    /* } */

    /* printf("\nThe elements of matrix for id %d are:\n",rank); */
    /* for(i=0;i<nrows;i++) */
    /* { */
    /* 	for(j=0;j<n;j++)  */
    /* 	{ */
    /* 	    printf("%0.0f  ",(*matrix)[i][j]); */
    /* 	}	        */
    /* 	printf("\n"); */
    /* } */

    /* MPI_Barrier(MPI_COMM_WORLD); */

    /* /\* Begin time measurement *\/ */
    /* if(rank==0) */
    /* { */
    /* 	gettimeofday(&startTime,NULL); */
    /* } */


    /* Do the matrix vector product in parallel */
    /* Each process simply performs an inner product */
    for(i=offset;i<nrows+offset;i++)
    {
	res[i] = 0.;
    	for(j=0;j<n;j++)
	{
    	    res[i]+=vector[j]*(*matrix)[i-offset][j];
	}
    }
    
    /* All processes send their result back to process 0,
       which puts the results in the result vector and prints them */
    if(rank==0)
    {
    	for(i=1;i<p;i++)
    	{
    	    MPI_Irecv(&res[displ[i]],length[i],MPI_DOUBLE,i,i,
    		     MPI_COMM_WORLD,&req[i-1]);
    	}
	MPI_Waitall(p-1,&req[0],MPI_STATUSES_IGNORE);

	/* Print result vector */
    	//printf("\n\nThe result vector is: \n\n");
    	//for(i=0;i<n;i++) { printf("  %f\n",res[i]); }

	/* Get time */
	gettimeofday(&endTime,NULL);
	double time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	printf("Matrix size: %d\t Number of processes: %d\t Time: %f s\n",n,p,time);
    }
    else
    {
    	MPI_Send(&res[offset],nrows,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
    }

    /* Terminate MPI */
    MPI_Finalize();
    
    /* Cleaning up */
    free(matrix);
    free(vector);
    free(res);
    return(0);
}
