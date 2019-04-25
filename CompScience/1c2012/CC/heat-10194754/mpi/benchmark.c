#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>

int main(int argc, char *argv[])
{
    //Initialization
    int num, id;

    MPI_Init(NULL,NULL);

    MPI_Comm_size(MPI_COMM_WORLD,&num);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);

    //# iterations and message size:
    //Latency: use large niter and small size (1)
    //Bandwith: use small niter and large size 
    int niter = atoi(argv[1]);
    int size = atoi(argv[2]);

    int i;
    int tag = 99;
    char *a = (char*) malloc(sizeof(char)*size);
    
    MPI_Status status;
    
    struct timeval startTime;
    struct timeval endTime;

    MPI_Request *req = (MPI_Request*) malloc(sizeof(MPI_Request)*2);

    //Message sending and receiving initialization
    if(id==0)
    {
	MPI_Send_init(a,size,MPI_CHAR,1,tag,MPI_COMM_WORLD,&req[0]);
	MPI_Recv_init(a,size,MPI_CHAR,num-1,tag,MPI_COMM_WORLD,&req[1]);
    }
    else if(id==(num-1))
    {
	MPI_Send_init(a,size,MPI_CHAR,0,tag,MPI_COMM_WORLD,&req[0]);
	MPI_Recv_init(a,size,MPI_CHAR,num-2,tag,MPI_COMM_WORLD,&req[1]);
    }
    else
    {
	MPI_Send_init(a,size,MPI_CHAR,id+1,tag,MPI_COMM_WORLD,&req[0]);
	MPI_Recv_init(a,size,MPI_CHAR,id-1,tag,MPI_COMM_WORLD,&req[1]);	
    }

    //Barrier to synchronise all processes
    MPI_Barrier(MPI_COMM_WORLD);

    //Start time measurement
    if(id==0)
    {
	gettimeofday(&startTime,NULL);
    }

    //Perform iterations for the main process
    if(id==0)
    {

	for(i=0; i<niter; i++)
	{
	    MPI_Start(&req[0]);
	    MPI_Wait(&req[0],&status);
	    MPI_Start(&req[1]);
	    MPI_Wait(&req[1],&status);
	}
    }

    //Perform iterations for all other processes
    else
    {
	for(i=0; i<niter; i++)
	{
	    MPI_Start(&req[1]);
	    MPI_Wait(&req[1],&status);
	    MPI_Start(&req[0]);
	    MPI_Wait(&req[0],&status);
	}
    }
    
    //Calculate elapsed time and latency/bandwith
    if(id==0)
    {
	gettimeofday(&endTime,NULL);
	double time = (endTime.tv_sec - startTime.tv_sec) + 
	    ( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000;
	printf("It = %d, time = %f\n",niter,time);
	printf("Latency: %e seconds\n",time/(niter*num));
	double b = (double)size*num*niter;
	printf("Bandwith: %e MB per sec\n",b/(1000000.*time));
	printf("NOTE: latency is only accurate when the data size is 1!\n");
	printf("NOTE: bandwidth is only accurate if the data size is larger then 10^7 bytes!\n");
    }

    MPI_Finalize();
    free(a);

    return(1);
}
