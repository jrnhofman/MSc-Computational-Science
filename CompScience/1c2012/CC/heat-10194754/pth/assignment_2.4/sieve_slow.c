#include <stdio.h>
#include <stdlib.h>
#include "pthread.h"
#include <sys/time.h>

int BUFSIZE;
int primes;

int *buffer;
int thread_id;
pthread_t *thread_ids;
pthread_attr_t attr;

/* Generator function */
void *generator(void*a)
{
    int integer = 2;
    int in = 0;
   
    while(thread_id<primes)
    {
	if(buffer[in]==0)
	{
	    buffer[in] = integer;
	    in = (in + 1)%BUFSIZE;
	    integer++;
	}
    }

    //printf("Terminating generator!\n");
    return(NULL);
}

/* Filter function */
void *filter(void *a)
{

    int id = thread_id;

    int prime = *(int*)a;
    int candidate;
    int out = (id)*BUFSIZE;
    int in = (id-1)*BUFSIZE;
    int spawned = 0;

    printf("Prime %d: %d\n",thread_id,prime);

    do
    {
	if(buffer[out]==0 && buffer[in]!=0)
	{
	    candidate = buffer[in];
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
		    buffer[out] = candidate;
		    out = (out+1)%BUFSIZE + (id)*BUFSIZE;
		}
	    }
	    
	    buffer[in] = 0;
	    in = (in+1)%BUFSIZE + (id-1)*BUFSIZE;

	}
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

    /* Buffer, pthread initialisation */
    buffer = (int*) malloc(sizeof(int)*primes*BUFSIZE);

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

    /* Pthread creation */
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
    return(0);
}
