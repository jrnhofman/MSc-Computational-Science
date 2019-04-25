#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#define PI 3.1415926535897932384626433832795028

FILE *fc;
FILE *fp;

//Computing the call option price
void sameseed(int n,double S0, double K, double r, double sigma, double t, double epsilon)
{
    double std,error,stdb,errorb,payoff;
    double nrml,factor1,factor2;
    const double term = S0 * exp((r-sigma*sigma/2)*(1-t));
    const double termb = (S0+epsilon) * exp((r-sigma*sigma/2)*(1-t));
    const double terms = sigma*sqrt(1-t);
    double delta;
    int i;

    /*
     * Unbumped
     */

    //Seed
    srand(34874988);

    double var = 0.;
    double mean = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	payoff = fmax(0.,term*exp(nrml*terms)-K);

	//Dynamic mean/std
	var += (i-1)*(payoff-mean)*(payoff-mean)/i;
	mean += (payoff-mean)/i;

	i++;

	//The other random normal variable
	nrml = factor1*sin(factor2);

	//Call payoff
	payoff = fmax(0.,term*exp(nrml*terms)-K);

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
     * Bumped
     */

    //Seed
    srand(34874988);

    double varb = 0.;
    double meanb = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	payoff = fmax(0.,termb*exp(nrml*terms)-K);

	//Dynamic mean/std
	varb += (i-1)*(payoff-meanb)*(payoff-meanb)/i;
	meanb += (payoff-meanb)/i;

	i++;

	//The other random variable
	nrml = factor1*sin(factor2);

	//Call payoff
	payoff = fmax(0.,termb*exp(nrml*terms)-K);

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

    delta = (meanb - mean)/epsilon; 
    
    //printf("S0: %f, K: %f, r: %f, sigma: %f, t: %f, eps: %f\n",S0,K,r,sigma,t,epsilon);
    printf("Same seed\t N: %d\t epsilon: %f\t Delta: %f\n",n,epsilon,delta);
    //fprintf(fc,"%f %f %f\n",K,mean,error);
    //fprintf(fc,"%f %f %f\n",sigma,mean,error);

    //printf("European call (%d simulations): Mean value: %f, std. error: %f\n",n,mean,error);
}

//Computing the call option price
void diffseed(int n,double S0, double K, double r, double sigma, double t, double epsilon)
{
    double std,error,stdb,errorb,payoff;
    double nrml,factor1,factor2;
    const double term = S0 * exp((r-sigma*sigma/2)*(1-t));
    const double termb = (S0+epsilon) * exp((r-sigma*sigma/2)*(1-t));
    const double terms = sigma*sqrt(1-t);
    double delta;
    int i;

    /*
     * Unbumped
     */

    //Seed
    srand(34874988);

    double var = 0.;
    double mean = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	payoff = fmax(0.,term*exp(nrml*terms)-K);

	//Dynamic mean/std
	var += (i-1)*(payoff-mean)*(payoff-mean)/i;
	mean += (payoff-mean)/i;

	i++;

	//The other random normal variable
	nrml = factor1*sin(factor2);

	//Call payoff
	payoff = fmax(0.,term*exp(nrml*terms)-K);

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
     * Bumped
     */

    //Seed
    srand(64538299);

    double varb = 0.;
    double meanb = 0.;

    for(i=1; i<n+1; i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	payoff = fmax(0.,termb*exp(nrml*terms)-K);

	//Dynamic mean/std
	varb += (i-1)*(payoff-meanb)*(payoff-meanb)/i;
	meanb += (payoff-meanb)/i;

	i++;

	//The other random variable
	nrml = factor1*sin(factor2);

	//Call payoff
	payoff = fmax(0.,termb*exp(nrml*terms)-K);

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

    delta = (meanb - mean)/epsilon; 
    
    //printf("S0: %f, K: %f, r: %f, sigma: %f, t: %f, eps: %f\n",S0,K,r,sigma,t,epsilon);
    printf("Different seed\t N: %d\t epsilon: %f\t Delta: %f\n",n,epsilon,delta);
    //fprintf(fc,"%f %f %f\n",K,mean,error);
    //fprintf(fc,"%f %f %f\n",sigma,mean,error);

    //printf("European call (%d simulations): Mean value: %f, std. error: %f\n",n,mean,error);
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

    int i;
    for(i=5; i<=8; i++)
    {
	sameseed(pow(10,i),100.,99,0.06,0.20,0.,0.5);
	diffseed(pow(10,i),100.,99,0.06,0.20,0.,0.5);
	sameseed(pow(10,i),100.,99,0.06,0.20,0.,0.1);
	diffseed(pow(10,i),100.,99,0.06,0.20,0.,0.1);
	sameseed(pow(10,i),100.,99,0.06,0.20,0.,0.01);
	diffseed(pow(10,i),100.,99,0.06,0.20,0.,0.01);
	sameseed(pow(10,i),100.,99,0.06,0.20,0.,0.001);
	diffseed(pow(10,i),100.,99,0.06,0.20,0.,0.001);
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
    
