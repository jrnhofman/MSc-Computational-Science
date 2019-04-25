#include <stdio.h>
#include <stdlib.h>
#include "pthread.h"
#include <semaphore.h>
#include "sys/time.h"

int BUFSIZE;
int primes;

int *buffer;
int thread_id;
pthread_t *thread_ids;
pthread_attr_t attr;
sem_t *occupy;
sem_t *empty;

/* Generator function */
void *generator(void*a)
{
    int integer = 2;
    int in = 0;
   
    while(thread_id<primes)
    {
	sem_wait(&empty[0]);
	buffer[in] = integer;
	sem_post(&occupy[0]);
	in = (in + 1)%BUFSIZE;
	integer++;
    }
    //printf("Terminating generator!\n");
    return(NULL);
}

/* Filter functions */
void *filter(void *a)
{
    int id = thread_id;
    if(id==primes)
    {
	int i;
	for(i=0; i<primes; i++)
	{
	    sem_post(&empty[i]);
	    sem_post(&occupy[i]);
	}
    }

    int prime = *(int*)a;
    int candidate;
    int bufferbegin = (id-1)*BUFSIZE;
    int bufferend = id*BUFSIZE;
    int out = (id)*BUFSIZE;
    int in = (id-1)*BUFSIZE;
    int spawned = 0;

    printf("Prime %d: %d\n",thread_id,prime);

    do
    {
	sem_wait(&occupy[id-1]);
	candidate = buffer[in];
	sem_post(&empty[id-1]);
	if(candidate%prime!=0)
	{
	    if(spawned==0)
	    {
		int *num = (int*) malloc(sizeof(int));
		*num = candidate;
		thread_id++;
		pthread_create(&thread_ids[id],&attr,&filter,num);
		spawned = 1;
	    }
	    else
	    {
		sem_wait(&empty[id]);
		buffer[out] = candidate;
		sem_post(&occupy[id]);
		out = (out+1)%BUFSIZE + bufferend;
	    }
	}
	in = (in+1)%BUFSIZE + bufferbegin;

    } while(thread_id<primes); 
    //printf("Terminating filter %d!\n",id);
    return(NULL);
}

int main(int argc, char *argv[])
{
    if(argc==3) { primes = atoi(argv[1]); BUFSIZE = atoi(argv[2]); }
    if(argc==1)
    {
	printf("Enter amount of prime numbers you wish to compute:\n");
	if(scanf("%d",&primes)==1);
	printf("Enter buffer size:\n");
	if(scanf("%d",&BUFSIZE)==1);
    }


    int i;
    thread_id = 1;
    
    /* Buffer, pthread and semaphore initialisation */
    buffer = (int*) malloc(sizeof(int)*primes*BUFSIZE);
    occupy = (sem_t*) malloc(sizeof(sem_t)*primes);
    empty = (sem_t*) malloc(sizeof(sem_t)*primes);

    for(i=0; i<primes; i++)
    {
	sem_init(&occupy[i],0,0);
	sem_init(&empty[i],0,BUFSIZE);
    }

    pthread_t thread_gen;
    pthread_t threads[primes];
    thread_ids = threads;

    for(i=0; i<(primes*BUFSIZE); i++)
    {
	buffer[i] = 0;
    }

    if(pthread_attr_init(&attr))
    {
	printf("Creation of attributes failed!");
	exit(1);
    }

    pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);

    /* Timer and pthread creation */
    struct timeval startTime;
    struct timeval endTime;

    gettimeofday(&startTime,NULL);

    if(pthread_create(&thread_gen,&attr,&generator,NULL)!=0)
    {
	printf("Creation of pthread failed!");
	exit(1);
    }

    int *num = (int*) malloc(sizeof(int));
    *num = 2;     
    if(pthread_create(&thread_ids[0],&attr,&filter,num)!=0)
    {
	printf("Creation of pthread failed!");
	exit(1);
    }

    /* Thread joining */
    pthread_join(thread_gen,NULL);
    for(i=0; i<primes; i++)
    {
	pthread_join(thread_ids[i],NULL);
    }

    /* Time reporting */
    gettimeofday(&endTime,NULL);
    double tS = startTime.tv_sec + (double) startTime.tv_usec/1000000;
    double tE = endTime.tv_sec + (double) endTime.tv_usec/1000000;
    double time = tE - tS;
    printf("Time = %f\n",time);

    /* Cleaning up */
    pthread_attr_destroy(&attr);
    free(buffer);
    free(occupy);
    free(empty);
    return(0);
}
