#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#define PI 3.1415926535897932384626433832795028

FILE *fc;
FILE *fp;

//Computing the asian option price
void computec(int n,double S0, double K, double r, double sigma)
{
    double std,error,payoff,payoffinv;
    double nrml,factor1,factor2;
    double var = 0.;
    double mean = 0.;
    double S,meanS,Sinv,meanSinv;
    int i,t;

    const double factor = exp((r-sigma*sigma/2)*(1./365));
    const double squareroot = sqrt(1./365);

    //Seed
    srand(time(NULL));

    for(i=1; i<(n/2+1); i++)
    {
	meanS = 0;
	meanSinv = 0;

	S = S0;
	Sinv = S0;

	for(t=1; t<366; t++)
	{
	    meanS += (S-meanS)/t;
	    meanSinv += (Sinv-meanSinv)/t;

	    //Box Muller normally distributed random variables
	    factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	    factor2 = 2*PI*(double)rand()/RAND_MAX;
	    nrml = squareroot*factor1*cos(factor2);

	    S = S*factor*exp(sigma*nrml);
	    Sinv = Sinv*factor*exp(sigma*nrml*-1);

	    t++;

	    meanS += (S-meanS)/t;
	    meanSinv += (Sinv-meanSinv)/t;

	    //Box Muller normally distributed random variables
	    nrml = squareroot*factor1*sin(factor2);

	    S = S*factor*exp(sigma*nrml);
	    Sinv = Sinv*factor*exp(sigma*nrml*-1);

	}

	//Call payoff
	payoff = fmax(0.,meanS-K);
	payoffinv = fmax(0.,meanSinv-K);

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
    //fprintf(fc,"%f %f %f\n",sigma,mean,error);

    printf("S0: %f, K: %f, r: %f, sigma: %f\n",S0,K,r,sigma);
    printf("Asian (%d simulations): Mean value: %f, std. error: %f\n",n,mean,error);
}

int main(int argc, char **argv)
{
    /* /\* Time measurement *\/ */
    /* struct timeval startTime; */
    /* struct timeval endTime; */
    /* gettimeofday(&startTime,NULL); */

    computec(100000,100.,99.,0.06,0.20);
    computec(1000000,100.,99.,0.06,0.20);
    computec(10000000,100.,99.,0.06,0.20);
    computec(100000000,100.,99.,0.06,0.20);

    /* fc = fopen("strike_call_asian.dat","w"); */

    /* int i; */

    /* for(i=0;i<201;i+=10) */
    /* { */
    /* 	computec(1000000,100.,(double)i,0.06,0.20); */
    /* } */

    /* fclose(fc); */
 

    /* fc = fopen("vol_call_asian.dat","w"); */

    /* int i; */

    /* for(i=0;i<101;i+=5) */
    /* { */
    /* 	computec(1000000,100.,99.,0.06,(double)i/100.); */
    /* } */

    /* fclose(fc); */
 
    /* /\* Time measurement *\/ */
    /* gettimeofday(&endTime,NULL); */
    /* double time = (endTime.tv_sec - startTime.tv_sec) +  */
    /* 	( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000; */
    /* printf("Execution Time: %f\n",time); */

    return(1);
}
    
