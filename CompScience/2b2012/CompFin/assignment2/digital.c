#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#define PI 3.1415926535897932384626433832795028

FILE *fc;
FILE *fp;
int seed;

//Computing the call option price
void euler(int n,double S0, double K, double r, double sigma, double epsilon)
{
    double std,error,stdb,errorb,stdc,errorc,payoff;
    double nrml,factor1,factor2;
    const double term = S0 * exp(r-sigma*sigma/2);
    const double termb = (S0+epsilon) * exp(r-sigma*sigma/2);
    const double termc = (S0-epsilon) * exp(r-sigma*sigma/2);
    const double terms = sigma;
    double delta,deltac;
    int i;

    seed++;

    /*
     * Unbumped
     */

    //Seed
    srand(seed);

    double var = 0.;
    double mean = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	if((term*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	//Dynamic mean/std
	var += (i-1)*(payoff-mean)*(payoff-mean)/i;
	mean += (payoff-mean)/i;

	i++;

	//The other random normal variable
	nrml = factor1*sin(factor2);

	//Call payoff
	if((term*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	//Dynamic mean/std
	var += (i-1)*(payoff-mean)*(payoff-mean)/i;
	mean += (payoff-mean)/i;
    }

    var /= n-1;
    std = sqrt(var);
    
    //Std error and mean
    error = 1.96*std/sqrt(n);
    mean = mean*exp(-r);

    //printf("Unbumped: Mean: %f, Error: %f\n",mean,error);

    /*
     * Bumped (S0 + epsilon)
     */

    //Seed
    srand(seed);

    double varb = 0.;
    double meanb = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	if((termb*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	//Dynamic mean/std
	varb += (i-1)*(payoff-meanb)*(payoff-meanb)/i;
	meanb += (payoff-meanb)/i;

	i++;

	//The other random variable
	nrml = factor1*sin(factor2);

	//Call payoff
	if((termb*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	//Dynamic mean/std
	varb += (i-1)*(payoff-meanb)*(payoff-meanb)/i;
	meanb += (payoff-meanb)/i;

    }

    varb /= n-1;
    stdb = sqrt(varb);
    
    //Std error and mean
    errorb = 1.96*stdb/sqrt(n);
    meanb = meanb*exp(-r);

    //printf("Bumped: Mean: %f, Error: %f\n",meanb,errorb);

    /*
     * Bumped (S0 - epsilon)
     */

    //Seed
    srand(seed);

    double varc = 0.;
    double meanc = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	if((termc*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	//Dynamic mean/std
	varc += (i-1)*(payoff-meanc)*(payoff-meanc)/i;
	meanc += (payoff-meanc)/i;

	i++;

	//The other random variable
	nrml = factor1*sin(factor2);

	//Call payoff

	if((termc*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	//Dynamic mean/std
	varc += (i-1)*(payoff-meanc)*(payoff-meanc)/i;
	meanc += (payoff-meanc)/i;

    }

    varc /= n-1;
    stdc = sqrt(varc);
    
    //Std error and mean
    errorc = 1.96*stdc/sqrt(n);
    meanc = meanc*exp(-r);

    //printf("Bumped: Mean: %f, Error: %f\n",meanc,errorc);
    //printf("diff: %f\n",meanb-meanc);
    delta = (meanb - mean)/epsilon; 
    deltac = (meanb - meanc)/(2*epsilon);

    //printf("S0: %f, K: %f, r: %f, sigma: %f, t: %f, eps: %f\n",S0,K,r,sigma,t,epsilon);
    printf("Forward Euler:\t N: %d\t epsilon: %f\t Delta: %f\t Difference: %f\n",n,epsilon,delta,100.*(delta- 0.018206369779490)/0.018206369779490);
    printf("Central Euler:\t N: %d\t epsilon: %f\t Delta: %f\t Difference: %f\n",n,epsilon,deltac,100.*(deltac- 0.018206369779490)/0.018206369779490);
    //fprintf(fc,"%f %f %f\n",K,mean,error);
    //fprintf(fc,"%f %f %f\n",sigma,mean,error);

}

//Computing the call option price
void likelihood(int n,double S0, double K, double r, double sigma)
{
    double std,error,payoff;
    double nrml,factor1,factor2;
    const double term = S0 * exp(r-sigma*sigma/2);
    const double terms = sigma;
    const double factor = exp(-r)/(sigma*S0);
    double delta;
    int i;

    seed++;

    srand(seed);

    double var = 0.;
    double mean = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	if((term*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	delta = factor*payoff*nrml;

	//Dynamic mean/std
	var += (i-1)*(delta-mean)*(delta-mean)/i;
	mean += (delta-mean)/i;

	i++;

	//The other random normal variable
	nrml = factor1*sin(factor2);

	//Call payoff
	if((term*exp(nrml*terms)) > K)
	{
	    payoff = 1.;
	}
	else
	{
	    payoff = 0.;
	}

	delta = factor*payoff*nrml;

	//Dynamic mean/std
	var += (i-1)*(delta-mean)*(delta-mean)/i;
	mean += (delta-mean)/i;
    }

    var /= n-1;
    std = sqrt(var);
    
    //Std error and mean
    error = 1.96*std/sqrt(n);
    mean = mean;

    //printf("Unbumped: Mean: %f, Error: %f\n",mean,error);
    printf("N: %d\t Delta estim.: %f\t error: %f\t difference: %f\n",n,mean,error,100.*(mean- 0.018206369779490)/0.018206369779490);

}

int main(int argc, char **argv)
{
    /* /\* Time measurement *\/ */
    /* struct timeval startTime; */
    /* struct timeval endTime; */
    /* gettimeofday(&startTime,NULL); */

    /* int i; */
    /* for(i=0;i<365;i+=5) */
    /* { */
    /* 	computec(1000000,100.,99,0.06,0.20,(double)i/365.,0.0001); */
       
    /* } */

    seed = time(NULL);

    int i;

    for(i=5; i<=8; i++)
    {
	likelihood(pow(10,i),100.,99,0.06,0.20);
    }

    for(i=5; i<=8; i++)
    {
    	euler(pow(10,i),100.,99,0.06,0.20,0.5);
    	euler(pow(10,i),100.,99,0.06,0.20,0.1);
    	euler(pow(10,i),100.,99,0.06,0.20,0.01);
    	euler(pow(10,i),100.,99,0.06,0.20,0.001);
    }

    /* fc = fopen("strike_call.dat","w"); */
    /* fp = fopen("strike_put.dat","w"); */

    /* int i; */

    /* for(i=0;i<201;i+=10) */
    /* { */
    /* 	computec(1000000,100.,(double)i,0.06,0.20); */
    /* 	computep(1000000,100.,(double)i,0.06,0.20); */
    /* } */

    /* fclose(fc); */
    /* fclose(fp); */


    /* fc = fopen("vol_call.dat","w"); */
    /* fp = fopen("vol_put.dat","w"); */

    /* int i; */

    /* for(i=0;i<101;i+=5) */
    /* { */
    /* 	computec(1000000,100.,99.,0.06,(double)i/100.); */
    /* 	computep(1000000,100.,99.,0.06,(double)i/100.); */
    /* } */

    /* fclose(fc); */
    /* fclose(fp); */

    /* /\* Time measurement *\/ */
    /* gettimeofday(&endTime,NULL); */
    /* double time = (endTime.tv_sec - startTime.tv_sec) +  */
    /* 	( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000; */
    /* printf("Execution Time: %f\n",time); */

    return(1);
}
    
