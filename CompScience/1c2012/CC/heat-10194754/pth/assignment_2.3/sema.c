#include <stdio.h>
#include "pthread.h"
#include "sys/time.h"
#include <stdlib.h>
#include <semaphore.h>

/* Parameters, nthreads and count are usually the same */
const int maxIter = 800000;
int nthreads;
int count;
	
/* Funtions for semaphore approach */
struct sem_barrier_t
{
    sem_t lock1;
    sem_t lock2;
    sem_t lock3;
    int add;
    int maxCount;
};

struct sem_barrier_t barrier_s;

void mythread_barrier_init(struct sem_barrier_t *barrier, void *a, int max)
{
    sem_init(&barrier->lock1,0,1);
    sem_init(&barrier->lock2,0,0);
    sem_init(&barrier->lock3,0,1);
    barrier->add = 0;
    barrier->maxCount = max;
}

void mythread_barrier_wait(struct sem_barrier_t *barrier)
{
    
    sem_wait(&barrier->lock1);
    barrier->add++;
    if(barrier->add==barrier->maxCount)
    {
	sem_wait(&barrier->lock3);
	sem_post(&barrier->lock2);
    }
    sem_post(&barrier->lock1);

    sem_wait(&barrier->lock2);
    sem_post(&barrier->lock2);

    sem_wait(&barrier->lock1);
    barrier->add--;
    if(barrier->add==0)
    {
	sem_wait(&barrier->lock2);
	sem_post(&barrier->lock3);
    }
    sem_post(&barrier->lock1);
    
    sem_wait(&barrier->lock3);
    sem_post(&barrier->lock3);
    
}

void sem_barrier_destroy(struct sem_barrier_t *barrier)
{
    sem_destroy(&barrier->lock1);
    sem_destroy(&barrier->lock2);
    sem_destroy(&barrier->lock3);
}
	
void *sem(void *id)
{
    int i;
    for(i=0; i<maxIter; i++)
    {
	mythread_barrier_wait(&barrier_s);
    }
    return(NULL);
}
  

int main(int argc, char *argv[])
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

    /* Barrier Initialisation */
    mythread_barrier_init(&barrier_s,NULL,count);

    /* Thread creation and joining */
    for(i=0; i<nthreads; i++)
    {
	if(pthread_create(&thread_ids[i],NULL,&sem,NULL)!=0)
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

    printf("Semaphores:\n");
    printf("Elapsed time %f\n",time);
    printf("Elapsed time per iteration %f\n",time/maxIter);
    printf("\n");

    /* Cleaning up */
    sem_barrier_destroy(&barrier_s);
    pthread_attr_destroy(&attr);
    return(1);
}
