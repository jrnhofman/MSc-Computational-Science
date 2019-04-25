#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include <sys/time.h>

#define sqrt 1.414213562373095048
#define dirw (sqrt/(4*(sqrt+1)))
#define diagw (1/(4*(sqrt+1)))
#define MAKEARRAY_2D(Type,Name,N,M)				\
    Type (*restrict Name)[N][M] =				\
	(Type (*restrict)[N][M]) malloc((N)*(M)*sizeof(Type))

void do_compute(const struct parameters* p, struct results *r)
{
    //Initialise variables
    size_t row,column, niter;
    double tS,tE;
    double diff;
    double maxdiff;
    
    const size_t N = p->N; //rows
    const size_t M = p->M; //columns
    double (*restrict old)[N+2][M+2];

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
    **  Initialise currIt and nextIt (where necessary)
    */

    double (*restrict c)[N+2][M+2] = malloc((N+2)*(M+2)*sizeof(double));

    /* Middle part */
    for(row=1; row<N+1; row++)
    {
	for(column=1; column<M+1; column++)
	{
	    (*currIt)[row][column] = p->tinit[(row-1)*M + column -1];
	    (*c)[row][column] = p->conductivity[(row-1)*M + column - 1];
	}
    }
    
    /* Top and bottom edges (FIXED) */
    for(column=1; column<M+1; column++)
    {
	(*currIt)[0][column] = (*currIt)[1][column];
	(*nextIt)[0][column] = (*currIt)[1][column];
	(*currIt)[N+1][column] = (*currIt)[N][column];
	(*nextIt)[N+1][column] = (*currIt)[N][column];
    }

    /* Edges (FIXED) */
    (*currIt)[0][0] = (*currIt)[1][M]; 
    (*currIt)[0][M+1] = (*currIt)[1][1]; 
    (*currIt)[N+1][0] = (*currIt)[N][M]; 
    (*currIt)[N+1][M+1] = (*currIt)[N][1]; 
    (*nextIt)[0][0] = (*currIt)[1][M]; 
    (*nextIt)[0][M+1] = (*currIt)[1][1]; 
    (*nextIt)[N+1][0] = (*currIt)[N][M]; 
    (*nextIt)[N+1][M+1] = (*currIt)[N][1];     

    /*
    ** Start of timer
    */

    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime, NULL);

    /*
    ** Execution of iteration loop
    ** Note: r->niter starts at 1 because of FLOPS calculation! (../src/output.c)
    */

    for(niter=1; niter<p->maxiter+1; niter++)
    {

	/*
	** Implement halo columns for continuous boundary conditions
	*/
	
	for(row=1; row<N+1; row++)
	{
	    (*currIt)[row][0] = (*currIt)[row][M];
	    (*currIt)[row][M+1] = (*currIt)[row][1];
	}

	//Set difference to zero
	maxdiff = 0.;

	//Iteration part
	for(row=1; row<N+1; row++)
	{
	    for(column=1; column<M+1; column++)
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
		diff = (diff>fabs((*currIt)[row][column]-(*nextIt)[row][column])) ? diff: fabs((*currIt)[row][column]-(*nextIt)[row][column]); 			  
	    }
	}

	    //Exit if treshold is reached
	    //Increment counter to sync with natural loop ending
	    if(maxdiff<p->threshold)
	    {
		niter++;
		break;
	    }

	    //Swap arrays
	    old = currIt;
	    currIt = nextIt;
	    nextIt = old;

	    /*
	    ** Results (NOTE: ARRAYS ARE SWAPPED)
	    */
 
	    if(niter%p->period==0)
	    {
		r->tmax = (*currIt)[1][1];
		r->tmin = (*currIt)[1][1];
		r->tavg = 0.;
		r->niter = niter;
		r->maxdiff = maxdiff;

		//Calculate mean, max, min and maxdiff
		for(row=1; row<N+1; row++)
		{
		    for(column=1; column<M+1; column++)
		    {
			r->tavg += (*currIt)[row][column];
			r->tmax = fmax(r->tmax, (*currIt)[row][column]);
			r->tmin = fmin(r->tmin, (*currIt)[row][column]);
		    }
		}

		r->tavg /= N*M;

		//Calculate elapsed time
		gettimeofday(&endTime,NULL);
		tS = startTime.tv_sec + (double) startTime.tv_usec/1000000;
		tE = endTime.tv_sec + (double) endTime.tv_usec/1000000;
		r->time = tE - tS;

		//Report results
		report_results(p,r);
	    }
	
	}	


	//Decrement iteration counter to calculate proper FLOPS
	niter--;

	/*
	** Final results (NOTE: ARRAYS ARE SWAPPED)
	*/
	r->tmax = (*currIt)[1][1];
	r->tmin = (*currIt)[1][1];
	r->tavg = 0.;
	r->niter = niter;
	r->maxdiff = maxdiff;

	//Calculate mean, max, min and maxdiff
	for(row=1; row<N+1; row++)
	{
	    for(column=1; column<M+1; column++)
	    {
		r->tavg += (*currIt)[row][column];
		r->tmax = fmax(r->tmax, (*currIt)[row][column]);
		r->tmin = fmin(r->tmin, (*currIt)[row][column]);
	    }
	}

	r->tavg /= N*M;

	//Calculate elapsed time
	gettimeofday(&endTime,NULL);
	tS = startTime.tv_sec + (double) startTime.tv_usec/1000000;
	tE = endTime.tv_sec + (double) endTime.tv_usec/1000000;
	r->time = tE - tS;
    
	//Free arrays
	free(currIt);
	free(nextIt);
    }
