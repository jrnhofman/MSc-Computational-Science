#include <stdio.h>
#include "pthread.h"
#include "sys/time.h"
#include <stdlib.h>
#include <semaphore.h>

/* Parameters, nthreads and count are usually the same */
const int maxIter = 800000;
int nthreads;
int count;

/* Functions for original POSIX */
pthread_barrier_t barrier;

void *posix(void *id)
{
    int i;
    for(i=0; i<maxIter; i++)
    {
	pthread_barrier_wait(&barrier);
    }
    return(NULL);
}
  

int main(void)
{
    printf("Enter number of threads:\n");
    if(scanf("%d",&nthreads)==1);
    count = nthreads;

    if(nthreads < count)
    {
	printf("Less threads than barrier counts, not possible!\n");
	exit(1);
    }
    
    /* Initialisation of variables and attributes */
    struct timeval startTime;
    struct timeval endTime;

    double time,tS,tE;
    int i;

    pthread_t thread_ids[nthreads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);

    gettimeofday(&startTime,NULL);

    /* Barrier initialisation */
    pthread_barrier_init(&barrier,NULL,count);

    /* Thread creation and joining */
    for(i=0; i<nthreads; i++)
    {
	if(pthread_create(&thread_ids[i],NULL,&posix,NULL)!=0)
	{
	    printf("Creation of pthread failed!");
	    exit(1);
	}
    }

    for(i=0; i<nthreads; i++)
    {
	pthread_join(thread_ids[i],NULL);
    }

    gettimeofday(&endTime,NULL);
    tS = startTime.tv_sec + (double) startTime.tv_usec/1000000;
    tE = endTime.tv_sec + (double) endTime.tv_usec/1000000;
    time = tE - tS;

    printf("Original POSIX functions:\n");
    printf("Elapsed time %f\n",time);
    printf("Elapsed time per iteration %f\n",time/maxIter);
    printf("\n");
    
    /* Cleaning up */
    pthread_barrier_destroy(&barrier);
    pthread_attr_destroy(&attr);
    return(1);
}
