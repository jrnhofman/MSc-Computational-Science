#include <stdio.h>
#include "pthread.h"
#include "sys/time.h"
#include <stdlib.h>
#include <semaphore.h>

/* Parameters, nthreads and count are usually the same */
const int maxIter = 800000;
int nthreads;
int count;

/* Functions for mutex locks/conditions */
struct barrier_t 
{
    int max;
    int number;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
}; 

struct barrier_t barrier_a;

void mythread_barrier_init(struct barrier_t *barrier,void *a,int max)
{
    barrier->max = max;
    barrier->number = 0;
    pthread_mutex_init(&barrier->mutex,NULL);
    pthread_cond_init(&barrier->cond,NULL);
}

void mythread_barrier_wait(struct barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    barrier->number++;
    if(barrier->number == barrier->max)
    {
	barrier->number = 0;
	pthread_cond_broadcast(&barrier->cond);
    }
    else
    {
	pthread_cond_wait(&barrier->cond,&barrier->mutex);
    }
    pthread_mutex_unlock(&barrier->mutex);
}

void mythread_barrier_destroy(struct barrier_t *barrier)
{
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
}

void *mutexcond(void *id)
{
    int i;
    for(i=0; i<maxIter; i++)
    {
	mythread_barrier_wait(&barrier_a);
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
    mythread_barrier_init(&barrier_a,NULL,count);

    /* Thread creation and joining */
    for(i=0; i<nthreads; i++)
    {
	if(pthread_create(&thread_ids[i],NULL,&mutexcond,NULL)!=0)
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

    printf("Mutex locks and conditions:\n");
    printf("Elapsed time %f\n",time);
    printf("Elapsed time per iteration %f\n",time/maxIter);
    printf("\n");

    /* Cleaning up */
    mythread_barrier_destroy(&barrier_a);
    pthread_attr_destroy(&attr);
    return(1);
}
