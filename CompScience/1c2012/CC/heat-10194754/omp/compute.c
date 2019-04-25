#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "compute.h"
#include <sys/time.h>
#include "omp.h"

#define sqrt 1.414213562373095048
#define dirw (sqrt/(4*(sqrt+1)))
#define diagw (1/(4*(sqrt+1)))
#define MAKEARRAY_2D(Type,Name,N,M)		\
    Type (*restrict Name)[N][M] = \
	(Type (*restrict)[N][M]) malloc((N)*(M)*sizeof(Type))

void do_compute(const struct parameters* p, struct results *r)
{
    //Initialise variables
    size_t row,column,niter;
    double tS,tE;
    double localmax,localmin,diff;
    double avg = 0.;
    double maxdiff = p->threshold+1;
    double maxdiffx = p->threshold+1;
    r->tmax = -INFINITY;
    r->tmin = INFINITY;
    r->maxdiff = 0.;
    const double threshold = p->threshold;
    
    omp_set_num_threads(p->nthreads);

    const size_t N = p->N; //rows
    const size_t M = p->M; //columns
    const size_t N1 = N+1;
    const size_t M1 = M+1;
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

    /*Columns */
    for(row=1; row<N+1; row++)
    {
	(*currIt)[row][0] = (*currIt)[row][M];
	(*currIt)[row][M+1] = (*currIt)[row][1];
    }

    /*
    ** Start of timer
    */

    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime, NULL);

#pragma omp parallel \
    private(niter,row,column,diff,localmax,localmin)
    {
	/*
	** Main loop
	*/

      for(niter=1; (niter<(p->maxiter+1)) && (maxdiffx>threshold); niter++)
	{
	    diff = 0.;

	    //Main computation
#pragma omp for schedule(static) nowait
	    for(row=1; row<N1; row++)
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
		    diff = (diff>fabs((*currIt)[row][column]-(*nextIt)[row][column])) ? diff: fabs((*currIt)[row][column]-(*nextIt)[row][column]); 
		}
	    }

	    //Writing the local variable to the global variable
	    //The bug is in this part of the code, maxdiff = diff
	    //slows the code down by 50%
	    if(diff>maxdiff)
	    {
#pragma omp critical
		{
	    
	    if(diff > maxdiff) { maxdiff = diff; }
		}
	    }

	    //Results 
	    if(niter%p->period==0)
	    {
		localmax = -INFINITY;
		localmin = INFINITY;

		//Calculate mean, local max, local min 
#pragma omp for reduction(+:avg) nowait schedule(static)		
		for(row=1; row<N1; row++)
		{
		    for(column=1; column<M1; column++)
		    {
			avg = avg + (*nextIt)[row][column];
			localmax = fmax(localmax, (*nextIt)[row][column]);
			localmin = fmin(localmin, (*nextIt)[row][column]);
		    }
		}
		
		//Calculate global max and min
#pragma omp critical
		{
		    r->tmax = fmax(r->tmax,localmax);
		    r->tmin = fmin(r->tmin,localmin);
		}

	    }

	    //Not strictly necessary on LISA
	    //but one cannot assume the same threads
	    //are assigned the same chunks between two
	    //static assignments
#pragma omp barrier

	    //Copy columns
#pragma omp for schedule(static)
	    for(row=1; row<N1; row++)
	    {
		(*nextIt)[row][0] = (*nextIt)[row][M];
		(*nextIt)[row][M1] = (*nextIt)[row][1];
	    }

#pragma omp single 
	    {

		//Report results
		if(niter%p->period==0)
		{
		    r->niter = niter;
		    r->tavg = avg/(N*M);
		    r->maxdiff = maxdiff;

		    gettimeofday(&endTime,NULL);
		    tS = startTime.tv_sec + (double) startTime.tv_usec/1000000;
		    tE = endTime.tv_sec + (double) endTime.tv_usec/1000000;
		    r->time = tE - tS;

		    report_results(p,r);
		}

		//Swap arrays
		old = currIt;
		currIt = nextIt;
		nextIt = old;
	
		//Initialise variables for next iteration
		r->tmax = -INFINITY;
		r->tmin = INFINITY;
		maxdiffx = maxdiff;
		maxdiff = 0.;
		avg = 0.;
	    }
	    
	}	

	/*
	** Final results (NOTE: Arrays are swapped)
	*/

	localmax = -INFINITY;
	localmin = INFINITY;

	//Calculate mean, local max, local min 
#pragma omp for schedule(static) reduction(+:avg)		
	for(row=1; row<N1; row++)
	{
	    for(column=1; column<M1; column++)
	    {
		avg += (*currIt)[row][column];
		localmax = fmax(localmax, (*currIt)[row][column]);
		localmin = fmin(localmin, (*currIt)[row][column]);
	    }
	}
		
	//Calculate global max and min
#pragma omp critical
	{
	    if(r->tmax < localmax) { r->tmax = localmax; }
	    if(r->tmin > localmin) { r->tmin = localmin; }
	}

#pragma omp single nowait
	{
	    //Final results reporting
	    r->niter = niter-1;
	    r->tavg = avg/(N*M);
	    r->maxdiff = maxdiffx;

	    gettimeofday(&endTime,NULL);
	    tS = startTime.tv_sec + (double) startTime.tv_usec/1000000;
	    tE = endTime.tv_sec + (double) endTime.tv_usec/1000000;
	    r->time = tE - tS;
	report_results(p,r);
	    //Free arrays
	    free(currIt);
	    free(nextIt);
	}
    }
}
