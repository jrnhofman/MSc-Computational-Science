#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include <sys/time.h>
#include <mpi.h>

#define sqrt 1.414213562373095048
#define dirw (sqrt/(4*(sqrt+1)))
#define diagw (1/(4*(sqrt+1)))
#define MAKEARRAY_2D(Type,Name,N,M)				\
    Type (*restrict Name)[N][M] =				\
	(Type (*restrict)[N][M]) malloc((N)*(M)*sizeof(Type))

void do_compute(const struct parameters* p, struct results *r)
{
    int num;
    int id;
    int tag = 99;

    MPI_Init(NULL,NULL);

    MPI_Comm_size(MPI_COMM_WORLD,&num);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);

    MPI_Status status;
    MPI_Request req1,req2,req3,req4;

    //Initialise variables
    size_t row,column,niter;
    double localdiff,localmax,localmin,localavg;
    r->maxdiff = p->threshold+1; 
   
    const size_t N = (id < ((p->N)%num)) ? ((p->N)/num + 1) : ((p->N)/num); //rows
    const size_t M = p->M; //columns
    const size_t N1 = N + 1;
    const size_t M1 = M + 1;
    const size_t M2 = M + 2;
    const size_t id_1 = id-1;
    const size_t id1 = id+1;
    const size_t num_1 = num -1;
    const size_t offset = (id <= ((p->N)%num)) ? (id*((p->N)/num) + id) : (id*((p->N)/num) + (p->N)%num);
    void *old;

    struct timeval startTime;
    struct timeval endTime;

    //Make two temperature arrays
    MAKEARRAY_2D(double,currIt,N+2,M+2);
    if(currIt == NULL)
    {
    	printf("Could not allocate memory, terminating!\n");
    	exit(1);
    }
    MAKEARRAY_2D(double,nextIt,N+2,M+2);
    if(nextIt == NULL)
    {
    	printf("Could not allocate memory, terminating!\n");
    	exit(1);
    }

    /*
    **  Initialise currIt and nextIt per process
    */

    double (*restrict c)[N+2][M+2] = malloc((N+2)*(M+2)*sizeof(double));

    //Only 1 process
    if(num==1)
    {
	for(row=1; row<N1; row++)
    	{
    	    for(column=1; column<M1; column++)
    	    {
    		(*currIt)[row][column] = p->tinit[(row+offset-1)*M + column -1];
    		(*c)[row][column] = p->conductivity[(row+offset-1)*M + column - 1];
    	    }
    	}

    	for(column=1; column<M1; column++)
    	{
    	    (*currIt)[0][column] = (*currIt)[1][column];
    	    (*nextIt)[0][column] = (*currIt)[1][column];
    	    (*currIt)[N+1][column] = (*currIt)[N][column];
    	    (*nextIt)[N+1][column] = (*currIt)[N][column];
    	}
    	for(row=1; row<N1; row++)
    	{
    	    (*currIt)[row][0] = (*currIt)[row][M];
    	    (*currIt)[row][M+1] = (*currIt)[row][1];
    	}

        (*currIt)[0][0] = (*currIt)[1][0];
    	(*nextIt)[0][0] = (*currIt)[1][0];
    	(*currIt)[0][M+1] = (*currIt)[1][M+1];
    	(*nextIt)[0][M+1] = (*currIt)[1][M+1];
 	(*currIt)[N+1][0] = (*currIt)[N][0];
    	(*nextIt)[N+1][0] = (*currIt)[N][0];
    	(*currIt)[N+1][M+1] = (*currIt)[N][M+1];
    	(*nextIt)[N+1][M+1] = (*currIt)[N][M+1];
    }

    //More than 1 process
    else
    {

    	for(row=1; row<N1; row++)
    	{
    	    for(column=1; column<M1; column++)
    	    {
    		(*currIt)[row][column] = p->tinit[(row+offset-1)*M + column -1];
    		(*c)[row][column] = p->conductivity[(row+offset-1)*M + column - 1];
    	    }
    	}
    	for(row=1; row<N1; row++)
    	{
    	    (*currIt)[row][0] = (*currIt)[row][M];
    	    (*currIt)[row][M+1] = (*currIt)[row][1];
    	}

	//Top process
    	if(id==0)
    	{
    	    for(column=1; column<M1; column++)
    	    {
    		(*currIt)[0][column] = (*currIt)[1][column];
    		(*nextIt)[0][column] = (*currIt)[1][column];

    		(*currIt)[N+1][column] = p->tinit[N*M+column-1];
    	    }
	    (*currIt)[N+1][0] = (*currIt)[N+1][M];
	    (*currIt)[N+1][M+1] = (*currIt)[N+1][1];
    	    (*currIt)[0][0] = (*currIt)[1][0];
    	    (*nextIt)[0][0] = (*currIt)[1][0];
    	    (*currIt)[0][M+1] = (*currIt)[1][M+1];
    	    (*nextIt)[0][M+1] = (*currIt)[1][M+1];
    	}

	//Bottom process
    	else if(id==num-1)
    	{
    	    for(column=1; column<M1; column++)
    	    {
    		(*currIt)[N+1][column] = (*currIt)[N][column];
    		(*nextIt)[N+1][column] = (*currIt)[N][column];

		(*currIt)[0][column] = p->tinit[(offset-1)*M + column - 1];
    	    }
    	    (*currIt)[N+1][0] = (*currIt)[N][0];
    	    (*currIt)[0][0] = (*currIt)[0][M];
	    (*currIt)[0][M+1] = (*currIt)[0][1];
	    (*nextIt)[N+1][0] = (*currIt)[N][0];
    	    (*currIt)[N+1][M+1] = (*currIt)[N][M+1];
    	    (*nextIt)[N+1][M+1] = (*currIt)[N][M+1];
    	}

	//All others
    	else
    	{
    	    for(column=1; column<M1; column++)
    	    {
    		(*currIt)[0][column] = p->tinit[(0+offset-1)*M + column -1];
    		(*currIt)[N+1][column] = p->tinit[(N+1+offset-1)*M + column -1];
    	    }

    	    (*currIt)[0][0] = (*currIt)[0][M];
    	    (*currIt)[0][M+1] = (*currIt)[0][1];

    	    (*currIt)[N+1][0] = (*currIt)[N+1][M];
    	    (*currIt)[N+1][M+1] = (*currIt)[N+1][1];
    	}
    }
   
    //Start timer
    if(id==0)
    {
	gettimeofday(&startTime, NULL);
    }

    //Main loop
    for(niter=1; (niter<(p->maxiter)+1) && (r->maxdiff>p->threshold); niter++)
    {
    	//Set difference to zero
    	localdiff = 0.;
	r->maxdiff = 0.;

	//Iterate over first and last row and copy corresponding halo column elements
	for(column=1; column<M1; column++)
	{
    	
	    (*nextIt)[1][column] = (*c)[1][column] * (*currIt)[1][column] +
		((*currIt)[1][column+1] +
		 (*currIt)[1][column-1] +
		 (*currIt)[2][column] +
		 (*currIt)[0][column]) * (1-(*c)[1][column]) * dirw +
		((*currIt)[0][column-1] +
		 (*currIt)[0][column+1] +
		 (*currIt)[2][column+1] +
		 (*currIt)[2][column-1]) * (1-(*c)[1][column]) * diagw;
	    localdiff = (localdiff>fabs((*currIt)[1][column]-(*nextIt)[1][column])) ? localdiff : fabs((*currIt)[1][column]-(*nextIt)[1][column]);
    	
	    (*nextIt)[N1-1][column] = (*c)[N1-1][column] * (*currIt)[N1-1][column] +
		((*currIt)[N1-1][column+1] +
		 (*currIt)[N1-1][column-1] +
		 (*currIt)[N1][column] +
		 (*currIt)[N1-2][column]) * (1-(*c)[N1-1][column]) * dirw +
		((*currIt)[N1-2][column-1] +
		 (*currIt)[N1-2][column+1] +
		 (*currIt)[N1][column+1] +
		 (*currIt)[N1][column-1]) * (1-(*c)[N1-1][column]) * diagw;
	    localdiff = (localdiff>fabs((*currIt)[N1-1][column]-(*nextIt)[N1-1][column])) ? localdiff : fabs((*currIt)[N1-1][column]-(*nextIt)[N1-1][column]);
	}
	(*nextIt)[1][0] = (*nextIt)[1][M];
	(*nextIt)[1][M1] = (*nextIt)[1][1];
	(*nextIt)[N1-1][0] = (*nextIt)[N1-1][M];
	(*nextIt)[N1-1][M1] = (*nextIt)[N1-1][1];
	
	//Initiate sending of rows 
	if(id!=0)
	{
	    MPI_Isend(&(*nextIt)[1][0],M2,MPI_DOUBLE,id_1,tag,MPI_COMM_WORLD,&req1);
	    MPI_Irecv(&(*nextIt)[0][0],M2,MPI_DOUBLE,id_1,tag,MPI_COMM_WORLD,&req3);
	}
	if(id!=num_1)
	{
	    MPI_Isend(&(*nextIt)[N][0],M2,MPI_DOUBLE,id1,tag,MPI_COMM_WORLD,&req2);
	    MPI_Irecv(&(*nextIt)[N1][0],M2,MPI_DOUBLE,id1,tag,MPI_COMM_WORLD,&req4);
	}
	
    	//Main computation
    	for(row=2; row<N; row++)
    	{
    	    for(column=1; column<M1; column++)
    	    {
    		(*nextIt)[row][column] = (*c)[row][column] * (*currIt)[row][column] +
    		    ((*currIt)[row][column+1] +
    		     (*currIt)[row][column-1] +
    		     (*currIt)[row+1][column] +
    		     (*currIt)[row-1][column]) * (1-(*c)[row][column]) * dirw +
    		    ((*currIt)[row-1][column-1] +
    		     (*currIt)[row-1][column+1] +
    		     (*currIt)[row+1][column+1] +
    		     (*currIt)[row+1][column-1]) * (1-(*c)[row][column]) * diagw;
    		localdiff = (localdiff>fabs((*currIt)[row][column]-(*nextIt)[row][column])) ? localdiff : fabs((*currIt)[row][column]-(*nextIt)[row][column]);
    	    }
    	}

        //Copy columns
    	for(row=2; row<N; row++)
    	{
    	    (*nextIt)[row][0] = (*nextIt)[row][M];
    	    (*nextIt)[row][M1] = (*nextIt)[row][1];
    	}

	//Reduce to global maxdiff value, this also works as a barrier
	MPI_Allreduce(&localdiff,&r->maxdiff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);	

	//Start sending and receiving and wait for completion
	if(id!=0)
	{
	    MPI_Wait(&req1,&status);
	    MPI_Wait(&req3,&status);
    	}

	if(id!=num_1)
	{
	    MPI_Wait(&req2,&status);
	    MPI_Wait(&req4,&status);
	}

	//Reporting
	if(niter%p->period==0)
	{
	    localmax = -INFINITY;
	    localmin = INFINITY;
	    localavg = 0.;

	    for(row=1; row<N1; row++)
	    {
		for(column=1; column<M1; column++)
		{
		    localavg += (*nextIt)[row][column];
		    localmax = fmax(localmax, (*nextIt)[row][column]);
		    localmin = fmin(localmin, (*nextIt)[row][column]);
		}
	    }

	    //Reductions of quantities
	    MPI_Reduce(&localavg,&r->tavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	    MPI_Reduce(&localmax,&r->tmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	    MPI_Reduce(&localmin,&r->tmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
	    MPI_Allreduce(&localdiff,&r->maxdiff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);	

	    if(id==0)
	    {
		r->niter = niter;
    		gettimeofday(&endTime,NULL);
		r->time = (endTime.tv_sec - startTime.tv_sec) + 
		    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
		r->tavg /= ((p->N) * (p->M));
		report_results(p,r);

		r->tavg = 0.;
		r->tmax = -INFINITY;
		r->tmin = INFINITY;
		r->maxdiff = 0.;
	    }
	}

	//Swap arrays
	old = currIt;
	currIt = nextIt;
	nextIt = old;

    }

    //Final results
    localmax = -INFINITY;
    localmin = INFINITY;
    localavg = 0.;

    for(row=1; row<N1; row++)
    {
	for(column=1; column<M1; column++)
	{
	    localavg += (*currIt)[row][column];
	    localmax = fmax(localmax, (*currIt)[row][column]);
	    localmin = fmin(localmin, (*currIt)[row][column]);
	}
    }

    //Reduction of quantities
    MPI_Reduce(&localavg,&r->tavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&localmax,&r->tmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&localmin,&r->tmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

    if(id==0)
    {
	r->niter = niter-1;
	gettimeofday(&endTime,NULL);
	r->time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	r->tavg /= (p->N) * (p->M);
	report_results(p,r);
    }
  
    //Cleaning up
    MPI_Finalize();
    free(currIt);
    free(nextIt);
    free(c);
}
