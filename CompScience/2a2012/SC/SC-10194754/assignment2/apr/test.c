
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>

#include "apr.h"

void test_fit_line()
{
	double x[] = {1, 2, 3};
	double y[] = {10, 20, 30};
	double a, b;
	int rank = -1;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	double res = fit_line(3, x, y, &a, &b);
	
	printf("debug (thread %i): fitline gives: a = %g, b = %g, res = %g\n", rank, a, b, res);

	assert((long)(a+0.5) == 10); // i know what the answer should be
	assert((long)(b+0.5) == 0);
	assert((long)(res+0.5) == 0);
}

void test_get_args(int argc, char* argv[])
{
	int rank = -1;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int ret = getparms(argc, argv); // use getparms and implicitly parseargs

	/* print all the parameter info that this thread has */

	printf("debug (thread %i): ret = %i\n", rank, ret);

	printf("debug (thread %i): parms.freq = %i\n", rank, parms.freq);
	printf("debug (thread %i): parms.ntotal = %i\n", rank, parms.ntotal);
	printf("debug (thread %i): parms.delta = %g\n", rank, parms.delta);
	printf("debug (thread %i): parms.stime = %g\n", rank, parms.stime);
	printf("debug (thread %i): parms.location = %g\n", rank, parms.location);
	printf("debug (thread %i): parms.periods = %i\n", rank, parms.periods);
}

void test_visualize()
{
	int rank = -1;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double mydata[] = {(rank*3.0), (rank*3.0)+1, (rank*3.0)+2};

	/* three data points, three time steps*/ 

	FILE* myfp = visualize_init(); // open my file

	assert(myfp != NULL); // make sure opening file succeeded

	visualize(myfp, 3, mydata, 1); // write to file
	mydata[0] += 1.0 ; // just to distinguish the data in time steps
	mydata[1] += 1.0 ;
	mydata[2] += 1.0 ;
	visualize(myfp, 3, mydata, 2); // write to file
	mydata[0] += 1.0 ; // just to distinguish the data in time steps
	mydata[1] += 1.0 ;
	mydata[2] += 1.0 ;
	visualize(myfp, 3, mydata, 3); // write to file

	visualize_combine();

	/* 
		Now there should be a file called string.data with the following,
		for three processes:

		0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0
		1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
		2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0
	*/
}

int main(int argc, char* argv[])
{
	int rank = -1;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	test_fit_line();

	test_get_args(argc, argv);

	test_visualize();

	MPI_Finalize();

	return 0;
} 
