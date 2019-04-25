#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>

int main(int argc, char *argv[])
{
    int size,i,j,n;
    int iter;
    int threads;
    
    if(argc==1)
    {
	printf("Enter grid size:\n:");
	if(scanf("%d",&size)==1);
	printf("Enter number of threads:\n");
	if(scanf("%d",&threads)==1);
	omp_set_num_threads(threads);
	printf("Enter iterations:\n");
	if(scanf("%d",&iter)==1);
    }
    else
    {
	omp_set_num_threads(atoi(argv[3]));
	iter = atoi(argv[2]);
	size = atoi(argv[1]);
    }

    //Matrix creation
    int **matrix;
    matrix = (int **) malloc(sizeof(int*)*size);
    if(matrix==NULL)
    {
	printf("Could not create matrix, aborting.\n");
	exit(1);
    }
    for(i=0; i<size; i++)
    {
	matrix[i] = (int *) malloc(sizeof(int)*(size-i));
	if(matrix[i]==NULL)
	{
	    printf("Could not create matrix, aborting.\n");
	    exit(1);
	}
    }

    //Matrix initialisation
    for(i=0; i<size; i++)
    {
	for(j=0; j<(size-i); j++)
	{
	    matrix[i][j] = i + (j + i);
	}
    } 


    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime,NULL);

    //Main loop
#pragma omp parallel private(i,j,n) 
{
    for(n=0; n<iter; n++)
    {
#pragma omp for schedule(guided,3) 
	for(i=0; i<size; i++)
	{
	    for(j=0; j<(size-i); j++)
	    {
		matrix[i][j]++;
	    }
	}
    }
}

    gettimeofday(&endTime,NULL);
    double tS = startTime.tv_sec + (double) startTime.tv_usec/1000000;
    double tE = endTime.tv_sec + (double) endTime.tv_usec/1000000;
    double time = tE - tS;
    printf("N = %d, \t time = %f\n",size,time);

    free(matrix);
    
    return(1);
}
