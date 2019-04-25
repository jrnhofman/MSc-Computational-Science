#include "MT/mt19937.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*Parameters*/
#define N_steps 10001
#define N_walks 1000
#define N_samples 1000
#define expect_val (N_walks/4)
#define expect_val_tot (expect_val*N_samples)

void sort(float *arr);
float maximum_D(float *arr);

int main()
{

    /*Random number generator from GSL*/
    const gsl_rng_type * T;
    //T = gsl_rng_randu; /*RandU*/
    T = gsl_rng_r250; /*R250*/
    gsl_rng * r = gsl_rng_alloc(T); 
    gsl_rng_set (r, 9287348);
    
    /*Seeds for other RNG's*/
    init_genrand(9287348); 
    srand48(9287348); 
    srand(9287348); 
    srandom(9287348);

    /*File with Quantis numbers*/
    FILE *f_read;
    f_read = fopen("quantis.txt","r");

    FILE *f;
    f = fopen("norm_data.dat","w"); /*File with norm data*/
 
    int i,j,k,x,y,quad,quadrant_total[4] = {0};
    float D,norm_sq,temp;
    float chisquared_total,sample_mean,sample_std;
    float chisquared[N_samples],mean[N_samples];
    long rnd;
 
    /*Compute chi squared N_samples times*/
    for(j=0; j<N_samples; j++)
    {
	if(j%(N_samples/10)==0) printf("%d %% complete\n",(j*100)/(N_samples));
	chisquared[j] = 0;
	mean[j] = 0;
	int quadrant[4] = {0};

	/*Compute N_walks random walks*/
	for(i=0; i<N_walks; i++)
	{
	    x = 0, y = 0, quad = 0;

	    for(k=0; k<N_steps; k++)
	    {
		/*Lines for different number generators*/
		//fscanf(f_read,"%ld",&rnd); /*Quantis*/
		//rnd = random(); 
		rnd = genrand_int32(); /*Mersenne Twister*/ 
		//rnd = lrand48(); 
		//rnd = rand(); 
		//rnd = gsl_rng_get (r); /*RNG from GSL library*/

		if(rnd & 1) { x++; }
		else { x--; }
	
		if(rnd & 2) { y++; }
		else { y--; }
	    }
	    
	    if(x<0) quad += 1;
	    if(y<0) quad += 2;
    
	    norm_sq = x*x + y*y; /*Squared norm*/
	    mean[j] += (norm_sq-mean[j])/(i+1); /*Mean squared norm*/
	    quadrant[quad]++;
	}

	fprintf(f,"%f\n",mean[j]);

	for(i=0; i<4; i++) {
	    temp = (float) quadrant[i]-expect_val;
	    chisquared[j] += temp * temp/expect_val; /*Chi2 for sample j*/
	    //printf("%f\n",chisquared[j]);

	    quadrant_total[i] += quadrant[i]; 
	}
    }

    /*Mean and std over all samples */
    for(i=0; i<N_samples; i++) { sample_mean += mean[i]/N_samples; }
    for(i=0; i<N_samples; i++) 
    { 
	temp = mean[i]-sample_mean;
	sample_std += temp * temp/(N_samples-1); 
    }
    printf("Mean norm squared = %f, std norm squared = %f\n",sample_mean,sqrt(sample_std));
    printf("Expected norm squared = %f\n",2.*N_steps);

    /*Sort chi2 array with insertion sort*/
    sort(chisquared);

    /*Compute D*/
    f = fopen("chi2_data.dat","w");
    for(i=0; i<N_samples; i++)
    {	
	fprintf(f,"%f %f\n",chisquared[i],(float)(i+1)/N_samples);
	temp = (float) (i+1)/N_samples - gsl_sf_gamma_inc_P(3./2,chisquared[i]/2); 
	if(temp>D) { D = temp; }
	temp = gsl_sf_gamma_inc_P(3./2,chisquared[i]/2) - (float) i/N_samples;
	if(temp>D) { D = temp; }
    }
    printf("D = %f\n",D);

    /*Compute total chi2*/
    for(i=0; i<4; i++) 
    { 
	temp = (float) quadrant_total[i]-expect_val_tot;
	chisquared_total += temp * temp/expect_val_tot; 
    }
    printf("chi-squared total = %f\n",chisquared_total);

    fclose(f);
    fclose(f_read);
    gsl_rng_free(r);
}

/*Insertion sort*/
void sort(float *arr)
{
    int i,j,run;
    float a;
    for(j=1; j<N_samples; j++)
    {
	a = arr[j];
	i = j-1;
	run = 1;
	while(run){ 
	    if(arr[i]>a)
	    {
		arr[i+1] = arr[i];
		i--;
		if(i<0) run = 0;
	    }
	    else run = 0;
	}
	arr[i+1] = a;
    }
}
