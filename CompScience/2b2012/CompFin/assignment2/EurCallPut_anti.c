#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#define PI 3.1415926535897932384626433832795028

FILE *fc;
FILE *fp;

//Computing the call option price
void computec(int n,double S0, double K, double r, double sigma)
{
    double std,error,payoff,payoffinv;
    double nrml,factor1,factor2;
    double var = 0.;
    double mean = 0.;
    double term = S0 * exp(r-sigma*sigma/2);
    int i;

    //Seed
    srand(time(NULL));

    for(i=1; i<(n/2+1); i++)
    {
	//Box Muller normally distributed random variables
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Call payoff
	payoff = fmax(0.,term*exp(sigma*nrml)-K);
	payoffinv = fmax(0.,term*exp(sigma*-1*nrml)-K);

	//Dynamic mean/std
	factor1 = (payoff+payoffinv)/2;
	var += (i-1)*(factor1-mean)*(factor1-mean)/i;
	mean += (factor1-mean)/i;
    }

    var /= (n/2-1);
    printf("Variance: %f\n",var);
    std = sqrt(var);
    
    //Std error and mean
    error = 1.96*std/sqrt(n);
    mean = mean*exp(-r);

    //fprintf(fc,"%f %f %f\n",K,mean,error);
    fprintf(fc,"%f %f %f\n",sigma,mean,error);

    printf("S0: %f, K: %f, r: %f, sigma: %f\n",S0,K,r,sigma);
    printf("European call (%d simulations): Mean value: %f, std. error: %f\n",n,mean,error);
}

//Compute European put option price
void computep(int n,double S0, double K, double r, double sigma)
{
    double std,error,payoff,payoffinv;
    double nrml,factor1,factor2;
    double var = 0.;
    double mean = 0.;
    double term = S0 * exp(r-sigma*sigma/2);
    int i;

    //Seed
    srand(time(NULL));

    for(i=1; i<(n/2+1); i++)
    {
	//Box Muller normally distributed random numbers
        factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	factor2 = 2*PI*(double)rand()/RAND_MAX;
	nrml = factor1*cos(factor2);

	//Put option payoff
	payoff = fmax(0.,K-term*exp(sigma*nrml));
	payoffinv = fmax(0.,K-term*exp(sigma*-1*nrml));

	//Dynamic mean/std
	factor1 = (payoff+payoffinv)/2;
	var += (i-1)*(factor1-mean)*(factor1-mean)/i;
	mean += (factor1-mean)/i;
    }

    var /= (n/2-1);
    printf("Variance: %f\n",var);
    std = sqrt(var);

    //Std error/mean
    error = 1.96*std/sqrt(n);
    mean = mean*exp(-r);

    //fprintf(fp,"%f %f %f\n",K,mean,error);
    fprintf(fp,"%f %f %f\n",sigma,mean,error);

    printf("S0: %f, K: %f, r: %f, sigma: %f\n",S0,K,r,sigma);
    printf("European put (%d simulations): Mean value: %f, std. error: %f\n",n,mean,error);
}

int main(int argc, char **argv)
{
    /* /\* Time measurement *\/ */
    /* struct timeval startTime; */
    /* struct timeval endTime; */
    /* gettimeofday(&startTime,NULL); */


    /* computec(100000,100.,99.,0.06,0.20); */
    /* computec(1000000,100.,99.,0.06,0.20); */
    /* computec(10000000,100.,99.,0.06,0.20); */
    /* computec(100000000,100.,99.,0.06,0.20); */
    /* computep(100000,100.,99.,0.06,0.20); */
    /* computep(1000000,100.,99.,0.06,0.20); */
    /* computep(10000000,100.,99.,0.06,0.20); */
    /* computep(100000000,100.,99.,0.06,0.20); */

    /* fc = fopen("strike_call_anti.dat","w"); */
    /* fp = fopen("strike_put_anti.dat","w"); */

    /* int i; */

    /* for(i=0;i<201;i+=10) */
    /* { */
    /* 	computec(1000000,100.,(double)i,0.06,0.20); */
    /* 	computep(1000000,100.,(double)i,0.06,0.20); */
    /* } */

    /* fclose(fc); */
    /* fclose(fp); */


    /* fc = fopen("vol_call_anti.dat","w"); */
    /* fp = fopen("vol_put_anti.dat","w"); */

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
    
