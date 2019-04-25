#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include <math.h>

int main(int argc, char* argv[])
{

    //Initialization
    int niter = atoi(argv[1]);
    int size = atoi(argv[2]);
    char *a = (char*)malloc(sizeof(char)*size);

    struct timeval startTime;
    struct timeval endTime;

    int num,id,i,j,round,tag=99;
    MPI_Init(NULL,NULL);

    MPI_Comm_size(MPI_COMM_WORLD,&num);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);

    MPI_Status *status = (MPI_Status*) malloc(sizeof(MPI_Status)*(num-1));
    MPI_Status status2;
    MPI_Request *req = (MPI_Request*) malloc(sizeof(MPI_Request)*(num-1));
    
    //Calculate the number of levels of the tree
    //for cascade broadcasting
    int num_t = num;
    int maxround = 1;
    if(((double)num_t/2)>1)
    {
	maxround ++;
	num_t /= 2;
    } 
    
    MPI_Barrier(MPI_COMM_WORLD);

    //Start timer
    if(id==0)
    {
    	gettimeofday(&startTime,NULL);
    }

    int idt;
    int round_shift;
    //Cascade broadcasting
    for(j=0; j<niter; j++)
    {
    	for(round=0; round<maxround; round++)
    	{
	    round_shift = 1 << round;
    	    if(id<round_shift)
    	    {
    		idt = id+round_shift;
    		if(idt<num)
    		{
    		    MPI_Send(a,size,MPI_CHAR,idt,tag,MPI_COMM_WORLD);
    		}
    	    }
    	    else if(id<num && id<(1 << round << 1))
    	    {
		MPI_Recv(a,size,MPI_CHAR,id-round_shift,tag,MPI_COMM_WORLD,&status2);
    	    }
    	    MPI_Barrier(MPI_COMM_WORLD);
    	}
    }
    
    //Report elapsed time per iteration
    if(id==0)
    {
    	gettimeofday(&endTime,NULL);
    	double time = (endTime.tv_sec - startTime.tv_sec) +
    	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
    	printf("Time cascade version = %e\n",time/niter);
    }

    //Barrier for all processes
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Start timer
    if(id==0)
    {
	gettimeofday(&startTime,NULL);
    }

    //Simple broadcasting, i.e. master thread sends message
    //to all other threads and waits for completion of sending
    for(j=0; j<niter; j++)
    {
	if(id==0)
	{
	    for(i=1; i<num; i++)
	    {
		MPI_Isend(a,size,MPI_CHAR,i,tag,MPI_COMM_WORLD,&req[i-1]);
	    }
	    MPI_Waitall(num-1,req,status);
	}
	else
	{
	    MPI_Recv(a,size,MPI_CHAR,0,tag,MPI_COMM_WORLD,&status2);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    }

    //Report elapsed time per iteration
    if(id==0)
    {
	gettimeofday(&endTime,NULL);
	double time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	printf("Time simple broadcasting = %e\n",time/niter);
    }

    //Barrier for new method
    MPI_Barrier(MPI_COMM_WORLD);

    //Start timer
    if(id==0)
    {
	gettimeofday(&startTime,NULL);
    }

    //Original BCast implementation
    for(j=0; j<niter; j++)
    {
	MPI_Bcast(a,size,MPI_CHAR,0,MPI_COMM_WORLD);
    }

    //Report elapsed time per iteration
    if(id==0)
    {
	gettimeofday(&endTime,NULL);
	double time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	printf("Time original broadcasting = %e\n",time/niter);
    }

    //Cleaning up
    MPI_Finalize();
    free(req);

    return(1);
}
