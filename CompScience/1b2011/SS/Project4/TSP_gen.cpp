#include "SFMT/SFMT.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
using namespace std;

/* Parameters */
const int nCities = 101; // Number of cities
const int nRoutes = 10000; // Population size
const int c = 10; // 1/c is percentage selected for offspring
const int maxGen = 1000; // Number of iterations best fit is allowed to be same before termination
const int nPerfect = 1; // Number of elitists
const double a = 1./nCities; // Probability of mutation
const bool mut = true; // Mutation on/off
const int nSim = 1; // Number of simulations

vector< vector<double> > cost(vector< vector<int> > &route, vector< vector<double> > &distM);

int main()
{
    /* Vector initialisations */
    vector< vector<int> > routes (nRoutes,vector<int> (nCities));
    vector< vector<double> > fitness (nRoutes, vector<double> (2));
    vector< vector<int> > newroutes (nRoutes,vector<int> (nCities));
    vector< vector<int> > neighbor (4,vector<int> (2,0) );

    int l, b, i, j, k, n, parent1, parent2, minN, currCity, parentf1, parentf2;
    bool occupy3;
    double best, xdiff, ydiff, number, min;
    double mean, old_mean = 0., std, old_std = 0., old_mean_n = 0., mean_n, minimum = pow(10,9);

    /* Save city coordinates */
    ofstream file;
    file.open("Data/pos_gen.dat");

    /* Save fitness as function of generation */
    ofstream file_data;
    file_data.open("Data/data_gen.dat");

    /* Save optimal route */
    ofstream file_out;
    file_out.open("Data/route_gen.dat");

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
	best = 0.;

	/* Initialize random routes (parents) */
	for(j=0; j<nRoutes; j++)
	{
	    for(i=0; i<nCities; i++) { routes[j][i] = i; }
	    for(i=0; i<nCities; i++) 
	    { 
		number = ((gen_rand32() % (nCities-i)) + i);
		swap(routes[j][i],routes[j][number]);
	    }
	}
	
	/* Calculate fitness of parents and sort */
	fitness = cost(routes,distM);
	sort(fitness.begin(),fitness.end());

	/* Genetic algorithm steps */
	n = 0;
	while(b<maxGen)
	{

	    /* Put best parent in offspring straight away (elitism) */
	    for(i=0; i<nCities; i++) { newroutes[0][i] = routes[fitness[0][1]][i]; }
	
	    /* Pick parents and compute offspring */
	    for(j=nPerfect; j<nRoutes; j++)
	    {
	    
		/* Pick parents random from c best ones */
		parentf1 = gen_rand32() % (nRoutes/c);
		parentf2 = gen_rand32() % (nRoutes/c);
		parent1 = fitness[parentf1][1];
		parent2 = fitness[parentf2][1];
		
		/* If parent fitness is the same, make a clone and skip Greedy algorithm */
		if(fitness[parentf1][0]==fitness[parentf2][0])
		{
		    for(i=0; i<nCities; i++) { newroutes[j][i] = routes[parent1][i]; }
		}
	       
		/* Else do Greedy algorithm */
		else
		{
		    /* Start with random city */
		    currCity = gen_rand32() % nCities;
		    newroutes[j][0] = currCity;

		    /* Create new route */
		    for(l=1; l<nCities; l++)
		    {

			/* Calculate neighbors of the current city */
			for(i=0; i<nCities; i++)
			{
			    if(routes[parent1][i]==currCity) { 
				neighbor[0][0] = routes[parent1][(i+1)%nCities]; 
				neighbor[1][0] = routes[parent1][(i-1)%nCities]; 
				break;
			    }
			}
			for(i=0; i<nCities; i++)
			{
			    if(routes[parent2][i]==currCity) { 
				neighbor[2][0] = routes[parent2][(i+1)%nCities]; 
				neighbor[3][0] = routes[parent2][(i-1)%nCities];
				break;
			    }
			}

			/* Check if neighbors are already used in the new route */
			for(i=0; i<l; i++)
			{
			    if(newroutes[j][i]==neighbor[0][0]) { neighbor[0][1] = 1; }
			    if(newroutes[j][i]==neighbor[1][0]) { neighbor[1][1] = 1; }
			    if(newroutes[j][i]==neighbor[2][0]) { neighbor[2][1] = 1; }
			    if(newroutes[j][i]==neighbor[3][0]) { neighbor[3][1] = 1; }
			}

			/* Compute closest neighbor */
			min = 100000.;
			for(i=0;i<4;i++)
			{
			    if(neighbor[i][1]!=1)
			    {
				if(distM[neighbor[i][0]][currCity] < min)
				{
				    min = distM[neighbor[i][0]][currCity];
				    minN = neighbor[i][0];
				}
			    }
			    else { neighbor[i][1] = 0; }
			}
		
			/* If not all neighbors are already used, add closest neighbor */
			if(min<99999.)
			{
			    newroutes[j][l] = minN;
			    currCity = minN;
			}
		
			/* Else choose a city randomly and add to new route */
			else 
			{
			    occupy3 = true;
			    while(occupy3)
			    {
				occupy3 = false;
				currCity = gen_rand32() % nCities;
				for(i=0; i<l; i++) { 
				    if(newroutes[j][i]==currCity) { occupy3 = true; break; } 
				}
			    }
			    newroutes[j][l] = currCity;
			}
		    }
		}

		/* Mutation */
		if(mut)
		{
		    for(i=0; i<nCities; i++)
		    {
			if(genrand_real1()<a)
			{
			    swap(newroutes[j][i],newroutes[j][gen_rand32() % nCities]);
			}
		    }
		}
	    }

	    /* After making new generation reset and compute fitness */
	    routes = newroutes;
	    fitness = cost(routes,distM);
	    sort(fitness.begin(),fitness.end());

            /* check if optimal solution is stuck */
	    if(fitness[0][0]==best) { b++; }
	    else
	    {
		b = 0;
		best = fitness[0][0];
	    }

	    file_data << fitness[0][0] << endl;
 
	    n++;
	}

	file_data.close();
	
	/* Data analysis */
	mean = old_mean + (fitness[0][0]-old_mean)/(k+1);
	if(k!=0) { std = (1-1./k)*old_std + (k+1.)*(old_mean-mean)*(old_mean-mean); }
	mean_n = old_mean_n + ((double)(n+1-maxGen)-old_mean_n)/(k+1);
	if(fitness[0][0] < minimum) { minimum = fitness[0][0]; }
	old_mean = mean;
	old_std = std;
	old_mean_n = mean_n;

	/* Print current best route with length */
	cout << "cost = " << fitness[0][0] << endl;
	cout << "route = " << endl;
	for(i=0; i<nCities; i++)
	{
	    cout << routes[fitness[0][1]][i] << "-"; 
	}
	cout << endl;
   
    }

    /* Print final statistics */
    cout << "Mean minima = " << mean << endl;
    cout << "Sample std  = " << sqrt(std) << endl;
    cout << "Avg steps = " << mean_n << endl;
    cout << "Min of minima = " << minimum << endl;
    
    /* Save optimal route */
    for(i=0; i<nCities; i++) { file_out << routes[fitness[0][1]][i] << endl; }
    file_out.close();
}
    
/* Calculate cost of routes */
vector< vector<double> > cost(vector< vector<int> > &routes, vector< vector<double> > &distM)
{
    int i, j;
    double sum;
    vector< vector<double> > fitness (nRoutes, vector<double> (2));
    for(j=0; j<nRoutes; j++)
    {
	sum = 0.;
	for(i=0; i<nCities-1; i++)
	{
	    sum += distM[routes[j][i]][routes[j][i+1]];
	}
	sum += distM[routes[j][nCities-1]][routes[j][0]];
	fitness[j][1] = j;
	fitness[j][0] = sum;
    } 
    return fitness;
}
