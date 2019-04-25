extern int foo(int,int);

int fib(int n)
{
    if(n<=1) return 1;
    else return fib(n-1) + fib(n-2);
}

void foobar(int k, int **b, int N, int M)
{
    int i,j,p=0,q;

#pragma omp parallel private(i,j) if(size>threshold)
    {

#pragma omp for schedule(static) 
	for(i=0; i<N; i++)
	{
	    for(j=0; j<M; j++)
	    {
		b[i][j] = i+j+k;
	    }
	}

#pragma omp for schedule(static) reduction(p:+)
	for(i=0; i<N; i++)
	{
	    for(j=0; j<M; j++)
	    {
		if(b[i][j]>k) { p = p+1; }
	    }
	}

#pragma omp single nowait
	printf("%d\n",p);

#pragma omp for schedule(guided)
	for(i=0; i<N; i++)
	{
	    for(j=0; j<M; j++)
	    {
		b[i][j] = fib(b[i][j]);
	    }
	}

    }

    for(i=1; i<N; i++)
    {
#pragma omp parallel for schedule(static) firstprivate(i) private(j) if(size>threshold)
	for(j=0; j<M; j++)
	{
	    b[i][j] = b[i][j] + b[i-1][j];
	}
    }
    

    for(i=0; i<N; i++)
    {
	for(j=0; j<M; j++)
	{
	    p = foo(b[i][j],p);
	}
    }

    printf("%d\n",p);
   
}
