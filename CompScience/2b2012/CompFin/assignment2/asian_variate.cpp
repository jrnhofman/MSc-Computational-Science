#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <boost/math/distributions/normal.hpp>


#define PI 3.1415926535897932384626433832795028

FILE *fc;
FILE *fp;

double geomprice(double S0, double K, double r, double sigma)
{
    using boost::math::normal;
    normal norm;
    double sigmaz = sigma * sqrt((2.*365 + 1)/(6.*(365+1)));
    double rho = ((r-0.5*sigma*sigma) + sigmaz*sigmaz)/2;
	
    double d1 = (log(S0/K) + (rho+0.5*sigmaz*sigmaz)) / (sigmaz);
    double d2 = d1 - sigmaz;
    double a = exp(-r)*(S0 * exp(rho)* cdf(norm,d1) - K*cdf(norm,d2));
    //printf("sigmaz %f, rho %f, d1 %f, d2 %f\n",sigmaz,rho,d1,d2);

    return a;
}

double cestimate(double n, double S0, double K, double r, double sigma)
{
    //Ex & Ey & Exy & Vary
    srand(time(NULL));
    int i, t;
    double mean = 0.;
    double meanxy = 0.;
    double meany = 0.;
    double vary = 0;
    double meanSgeom,meanS,S;
    double factor1,factor2,nrml,payoff,payoffgeom;

    const double factor = exp((r-sigma*sigma/2)*(1./365));
    const double squareroot = sqrt(1./365);

    for(i=1; i<(n/10+1); i++)
    {
	meanSgeom = 0;
	meanS = 0;
	S = S0;
	
	for(t=1; t<366; t++)
	{
	    meanS += (S-meanS)/t;
	    meanSgeom += log(S);

	    //Box Muller normally distributed random variables
	    factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
	    factor2 = 2*PI*(double)rand()/RAND_MAX;
	    nrml = squareroot*factor1*cos(factor2);

	    S = S*factor*exp(sigma*nrml);

	    t++;

	    meanS += (S-meanS)/t;
	    meanSgeom += log(S);

	    //Box Muller normally distributed random variables
	    nrml = squareroot*factor1*sin(factor2);

	    S = S*factor*exp(sigma*nrml);
	}

	meanSgeom = exp(meanSgeom/366);

	//Call payoff
	payoff = fmax(0.,meanS-K);
	payoffgeom = fmax(0.,meanSgeom-K);

	//Ey variance
	vary += (i-1)*(payoffgeom-meany)*(payoffgeom-meany)/i;
	meany += (payoffgeom-meany)/i;

	//Dynamic mean/std
	mean += (payoff-mean)/i;
	meanxy += (payoff*payoffgeom-meanxy)/i;
    }
    double ex = mean;
    double ey = meany;
    double exy = meanxy;
    vary /= (n/10)-1;

    double cov = exy - ex*ey;
    double c = -cov/vary;
    printf("Estimate of c: %f\n",c);
    printf("ex: %f, exy: %f, vary: %f, cov: %f, c: %f\n",ex,exy,vary,cov,c);

    return c;
} 


//Computing the asian option price with control variates
void computec(int n,double S0, double K, double r, double sigma,double c)
{
    double std,error,payoff,payoffgeom;
    double nrml,factor1,factor2;
    double var = 0.;
    double mean = 0.;
    double S,meanS,meanSgeom;
    int i,t;

    double mean2 = 0.;

    const double factor = exp((r-sigma*sigma/2)*(1./365));
    const double squareroot = sqrt(1./365);

    const double geomestimate = geomprice(S0,K,r,sigma); 
    printf("geompriceestimate: %f\n",geomestimate);

    // //Seed
    // srand(time(NULL));

    // for(i=1; i<(n+1); i++)
    // {
    // 	meanS = 0;
    // 	meanSgeom = 0;

    // 	S = S0;
	
    // 	for(t=1; t<366; t++)
    // 	{
    // 	    meanS += (S-meanS)/t;
    // 	    meanSgeom += log(S);

    // 	    //Box Muller normally distributed random variables
    // 	    factor1 = sqrt(-2*log((double)rand()/RAND_MAX));
    // 	    factor2 = 2*PI*(double)rand()/RAND_MAX;
    // 	    nrml = squareroot*factor1*cos(factor2);

    // 	    S = S*factor*exp(sigma*nrml);

    // 	    t++;

    // 	    meanS += (S-meanS)/t;
    // 	    meanSgeom += log(S);

    // 	    //printf("t: %d, meanS: %f, Sgeom: %f\n",t,meanS,exp(meanSgeom/t));

    // 	    //Box Muller normally distributed random variables
    // 	    nrml = squareroot*factor1*sin(factor2);

    // 	    S = S*factor*exp(sigma*nrml);

    // 	}

    // 	meanSgeom = exp(meanSgeom/366);

    // 	//printf("meanS: %f, meanSgeom: %f\n",meanS,meanSgeom);

    // 	//Call payoff
    // 	payoff = fmax(0.,meanS-K);
    // 	payoffgeom = fmax(0.,meanSgeom-K);

    // 	payoff = payoff + c*(payoffgeom - exp(r)*geomestimate); 

    // 	//Dynamic mean/std
    // 	var += (i-1)*(payoff-mean)*(payoff-mean)/i;
    // 	mean += (payoff-mean)/i;

    // 	mean2 += (payoffgeom-mean2)/i;
    // }

    // //printf("mean geom est.: %f\n",exp(-r)*mean2);

    // var /= (n-1);
    // //printf("Variance: %f\n",var);
    // std = sqrt(var);
    
    // //Std error and mean
    // error = 1.96*std/sqrt(n);
    // mean = mean*exp(-r);

    // //fprintf(fc,"%f %f %f\n",K,mean,error);
    // fprintf(fc,"%f %f %f\n",sigma,mean,error);

    // printf("S0: %f, K: %f, r: %f, sigma: %f\n",S0,K,r,sigma);
    // printf("Asian (%e simulations): Mean value: %f, std. error: %f\n",(double)n,mean,error);
}

int main(int argc, char **argv)
{
    /* /\* Time measurement *\/ */
    /* struct timeval startTime; */
    /* struct timeval endTime; */
    /* gettimeofday(&startTime,NULL); */

    //test(100000,100.,99.,0.06,0.20);
    double c = cestimate(1000000,100.,99.,0.06,0.20);
    
    computec(100000,100.,99.,0.06,0.20,c);
    // computec(1000000,100.,99.,0.06,0.20,c);
    //computec(10000000,100.,99.,0.06,0.20,c);
    // computec(100000000,100.,99.,0.06,0.20,c);

    // fc = fopen("strike_call_asian_variate.dat","w");

    // int i;

    // for(i=0;i<201;i+=10)
    // {
    // 	computec(1000000,100.,(double)i,0.06,0.20,c);
    // }

    // fclose(fc);
 

    // fc = fopen("vol_call_asian_variate.dat","w");

    // int i;

    // for(i=0;i<101;i+=5)
    // {
    // 	computec(100000,100.,99.,0.06,(double)i/100,c);
    // }

    // fclose(fc);
 
    /* /\* Time measurement *\/ */
    /* gettimeofday(&endTime,NULL); */
    /* double time = (endTime.tv_sec - startTime.tv_sec) +  */
    /* 	( (double)endTime.tv_usec - (double)startTime.tv_usec )/1000000; */
    /* printf("Execution Time: %f\n",time); */

    return(1);
}
    
