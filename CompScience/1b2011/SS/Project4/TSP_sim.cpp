#include "SFMT/SFMT.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

/* Parameters */
const int nCities = 101; // Number of cities
const int nSteps = 100000; // Number of steps per Markov chain
const int tempMeas = 1000000; // Number of random routes for temperature estimate 
const int nSim = 1; // Number of simulations
const int maxGen = 500; // Number of generations before halting simulation

double cost(vector<int> &route, vector< vector<double> > &distM);

int main()
{
    /* Vector initialisations */
    vector<int> bestRoute (nCities);
    vector<int> oldRoute (nCities);
    vector<int> newRoute (nCities);
    
    bool lowered;
    int i, j, k, n, b;
    double xdiff, ydiff, accept, random1, random2, delta;
    double oldCost, newCost, number, maxT, T, best;
    double mean, old_mean = 0., std, old_std = 0., old_mean_n = 0., mean_n;

    /* Save city coordinates */
    ofstream file;
    file.open("Data/pos_sim.dat");

    /* Save cost as function of Markov chain */
    ofstream file_data;
    file_data.open("Data/data_sim.dat");

    /* Save optimal route */
    ofstream file_out;
    file_out.open("Data/route_sim.dat");

    /* Initialize random cities and make distance matrix */
    /*
    init_gen_rand(0);
    vector< vector<double> > cityPos (nCities,vector<double> (2));
    for(i=0; i<nCities; i++)
    {
	cityPos[i][0] = genrand_real1();
	cityPos[i][1] = genrand_real1();
	file << cityPos[i][0] << " " << cityPos[i][1] << endl;
    }
    file.close();

    vector< vector<double> > distM (nCities,vector<double> (nCities,0));
    for(i=0; i<nCities; i++)
    {
	for(j=0; j<nCities; j++)
	{
	    xdiff = cityPos[i][0]-cityPos[j][0];
	    ydiff = cityPos[i][1]-cityPos[j][1];
	    distM[i][j] = sqrt(xdiff*xdiff + ydiff*ydiff);
	}
    }
    */

    /* Read in from coordinates (.tsp files) */
    ifstream file_pos;
    file_pos.open("Data/eil101.tsp");
    vector< vector<double> > cityPos (nCities,vector<double> (2));
    for(i=0; i<nCities; i++)
    {
	file_pos >> j >> cityPos[i][0] >> cityPos[i][1];
	cout << j << " " << cityPos[i][0] << " " << cityPos[i][1] <<endl;
	file << cityPos[i][0] << " " << cityPos[i][1] << endl;
    }
    file_pos.close();
    file.close();

    vector< vector<double> > distM (nCities,vector<double> (nCities,0));
    for(i=0; i<nCities; i++)
    {
	for(j=0; j<nCities; j++)
	{
	    xdiff = cityPos[i][0]-cityPos[j][0];
	    ydiff = cityPos[i][1]-cityPos[j][1];
	    distM[i][j] = sqrt(xdiff*xdiff + ydiff*ydiff);
	}
    }


    /* Read in from distance matrix (.txt files) */    
    /*
    ifstream file_pos;
    file_pos.open("Data/dantzig42_d.txt");
    vector< vector<double> > distM (nCities,vector<double> (nCities,0));
    for(i=0; i<nCities; i++)
    {
	for(j=0; j<nCities; j++)
	{
	    file_pos >> output;
	    distM[i][j] = output;
	    cout << output << endl;
	}
    }
    file_pos.close();
    */

    /* Loop over simulations */
    init_gen_rand(time(NULL));
    for(k=0; k<nSim; k++)
    {
	cout << "Simulation " << k << endl;
	b = 0;
	best = 1000000;

        /* Initialize random route */
	for(i=0; i<nCities; i++) { oldRoute[i] = i; }
	for(i=0; i<nCities; i++) 
	{ 
	    number = ((gen_rand32() % (nCities-i)) + i);
	    swap(oldRoute[i],oldRoute[number]);
	}

	/* Calculate cost */
	oldCost = cost(oldRoute,distM);

	/* Determining start temperature */
	maxT = 0.;
	for(j=0; j<tempMeas; j++)
	{
	    for(i=0; i<nCities; i++) 
	    { 
		swap(oldRoute[i],oldRoute[((gen_rand32() % (nCities-i)) + i)]);
	    }
	    newCost = cost(oldRoute,distM);
	    delta = abs(oldCost-newCost);
	    if(delta>maxT) { maxT = delta; }
       	    oldCost = newCost;
	}

	T = 2*maxT;
	cout << "Initial temperature = " << T << endl;
	
	/* Lower temperature until stop criterium */
	n = 0;
	while(b<maxGen)
	{
	    T *= 0.975;

	    /* Generate Markov chain */
	    lowered = false;
	    for(i=0; i<nSteps; i++)
	    {
		
		/* generate and apply lin2-opt transition */
		newRoute = oldRoute;
		random1 = gen_rand32() % nCities;
		random2 = gen_rand32() % nCities;
		if(random1<=random2) 
		{ 
		    for(j=0; j<=floor((random2-random1)/2); j++)
		    {
			swap(newRoute[random1+j],newRoute[random2-j]);
		    } 
		}
		else 
		{ 
		    for(j=0; j<=floor((random1-random2)/2); j++)
		    {
			swap(newRoute[random1-j],newRoute[random2+j]);
		    } 
		}

		/* Calculate new cost function */
		newCost = cost(newRoute,distM);

		/* Do Metropolis */
		if(newCost<oldCost) { oldRoute = newRoute; oldCost = newCost; }
		else 
		{
		    delta = newCost-oldCost;
		    accept = exp(-delta/T);
		    if(genrand_real1()<accept) { oldRoute = newRoute; oldCost = newCost; }
		}

		/* If a better route is found, save it */
		if(oldCost < best && i > nSteps/2) { 
		    best = oldCost;
		    bestRoute = oldRoute;
		    lowered = true;
		}

	    }

	    /* Write data every Markov chain */
	    file_data << oldCost << endl;
	    
	    /* check if optimal solution is stuck */
	    if(lowered) { b = 0; }
	    else { b++; }

	    n++;
	}

	file_data.close();

	/* Data analysis */
	mean = old_mean + (best-old_mean)/(k+1);
	if(k!=0) { std = (1-1./k)*old_std + (k+1.)*(old_mean-mean)*(old_mean-mean); }
	mean_n = old_mean_n + ((double)(n+1-maxGen)-old_mean_n)/(k+1);
	old_mean = mean;
	old_std = std;
	old_mean_n = mean_n;

	cout << "cost = " << best << endl;
	cout << "route = " << endl;
	for(i=0; i<nCities; i++)
	{
	    cout << bestRoute[i] << endl;
	}
	cout << endl;
	
    }

    /* Print final statistics */
    cout << "Mean minima = " << mean << endl;
    cout << "Std = " << std << endl;
    cout << "Error = " << sqrt(std) << endl;
    cout << "Avg steps = " << mean_n << endl;

    /* Save optimal route */
    for(i=0; i<nCities; i++) { file_out << bestRoute[i] << endl; }
    file_out.close();
}
 
/* Calculate cost of route */
    double cost(vector<int> &route, vector< vector<double> > &distM)
    {
	int i;
	double sum = 0.;
	for(i=0; i<nCities-1; i++)
	{
	    sum += distM[route[i]][route[i+1]];
	}
	sum += distM[route[nCities-1]][route[0]];

	return sum;
    }
