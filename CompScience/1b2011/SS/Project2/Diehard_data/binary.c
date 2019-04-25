#include "MT/mt19937.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    //T = gsl_rng_randu;
    T = gsl_rng_r250;
    r = gsl_rng_alloc (T);
    gsl_rng_set (r, 9287348);
    
    init_genrand(9287348); 
    srand48(9287348); 
    srand(9287348); 
    srandom(9287348);


    ifstream f_read;
    f_read.open("quantis.txt");
 
    int i;

    ofstream f_bin;
    f_bin.open("rand_data.dat");
    for(i=0;i<10000000;i++)
    {
	//long rnd;
	//f_read >> rnd;
	if(i%10==0) f_bin << "\n";
	//long rnd = genrand_int32();
	//long rnd = lrandom();
	//long rnd = lrand48(); 
	int rnd = rand(); 
	//long rnd = gsl_rng_get (r); 
	//cout << rnd << endl;
	f_bin << hex << abs(rnd);
    }
    f_bin.close();
}
