#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "MT/mt19937.h" 
using namespace std;

//Choose methods
const bool RS = false; //random sampling
const bool LH = true; //latin hypercube sampling
const bool OR = false; //orthogonal sampling

//Size and area of domain of sampling
const float MinRe = -2.0;
const float MaxRe = 1.0;
const float MinIm = -1.2;
const float MaxIm = 1.2;
const float Surface = (-MinRe+MaxRe)*(-MinIm+MaxIm);

//Constants
const int MaxIterations = 5000; //Maximum number of iterations
const int N = 10000; //Number of samples NOTE: MUST BE PERFECT SQUARE
const int N_sqrt = int(sqrt(float(N))); //Square root of number of samples
const int N_runs = 100; //Number of runs

//Scaling factor for hypercube and orthogonal sampling
const float x_Scale = (-MinRe + MaxRe)/N;
const float y_Scale = (-MinIm + MaxIm)/N;

//Function for random sampling
void random(vector<float> &Area_RS_It, vector<float> &Area_RS_Mean, vector<float> &Area_RS_Sq_Mean)
{
    int i,k;
    float x,y,Z_re,Z_im,Z_re2,Z_im2;
    
    //Loop over samples
    for(k=0; k<N; k++)
    {
	//Generate random samples
	x = genrand_real1()*(-MinRe+MaxRe)+MinRe;
	y = genrand_real1()*(-MinIm+MaxIm)+MinIm;
	Z_re = x, Z_im = y;
	
	//Loop over iterations
	for(i=0; i<=MaxIterations; i++)
	{
	    Z_re2 = Z_re*Z_re, Z_im2 = Z_im*Z_im;
	    if((Z_re2 + Z_im2) > 4) { break; }
	    Area_RS_It[i] += 1; //count samples that are still in Mandelbrot set at iteration i
	    Z_im = 2*Z_re*Z_im + y;
	    Z_re = Z_re2 - Z_im2 + x;
	}
    }
    //Calculate area as function of iteration and add to mean (squared) area
    for(i=0; i<(MaxIterations+1); i++) 
    { 
	Area_RS_It[i] *= (Surface/N); 
	Area_RS_Mean[i] += Area_RS_It[i];
	Area_RS_Sq_Mean[i] += Area_RS_It[i]*Area_RS_It[i];
    }
}

//Function for latin hypercube sampling
void latinhypercube(vector<float> &Area_LH_It, vector<float> &Area_LH_Mean, vector<float> &Area_LH_Sq_Mean)
{
    vector<int> Random_List (N);
    
    int i,k,number,old;
    float x,y,Z_re,Z_im,Z_re2,Z_im2;

    //Create random list, index = column number, value at index = row number
    for(i=0; i<N; i++) { Random_List[i] =  i; }
    for(i=0; i<N; i++) 
    { 
	number = (((long) genrand_int31() % (N-i)) + i), old = Random_List[i];
	Random_List[i] = Random_List[number];
	Random_List[number] =  old;
    }

    //Loop over samples
    for(k=0; k<N; k++)
    {
	//Generate random numbers according to random list
	x = (k + genrand_real2())*x_Scale + MinRe; //take columns left to right
	y = (Random_List[k] + genrand_real2())*y_Scale + MinIm; //take row according to random list
	Z_re = x, Z_im = y;

	//Loop over iterations
	for(i=0; i<=MaxIterations; i++)
	{
	    Z_re2 = Z_re*Z_re, Z_im2 = Z_im*Z_im;
	    if((Z_re2 + Z_im2) > 4) { break; }
	    Area_LH_It[i] += 1; //count samples that are still in Mandelbrot set at iteration i
	    Z_im = 2*Z_re*Z_im + y;
	    Z_re = Z_re2 - Z_im2 + x;
	}
    }
    //Calculate area as function of iteration and add to mean (squared) area
    for(i=0; i<(MaxIterations+1); i++) 
    { 
	Area_LH_It[i] *= (Surface/N); 
	Area_LH_Mean[i] += Area_LH_It[i];
	Area_LH_Sq_Mean[i] += Area_LH_It[i]*Area_LH_It[i];
    }
}

//Function for orthogonal sampling
void orthogonal(vector<float> &Area_OR_It, vector<float> &Area_OR_Mean, vector<float> &Area_OR_Sq_Mean)
{
    vector< vector<int> > xlist (N_sqrt,vector<int>(N_sqrt));
    vector< vector<int> > ylist (N_sqrt,vector<int>(N_sqrt));
    int i,j,number,old,n;
    int m = 0;
    float x,y,Z_re,Z_im,Z_re2,Z_im2;

    //Dividing area in N subsquares with each N rows and columns
    //First index of xlist,ylist corresponds to subsquare, second to column/row
    for(i=0; i<N_sqrt; i++)
    {
        for (j=0; j<N_sqrt; j++) { xlist[i][j] = ylist[i][j] = m++; }
    }
    
    //xlist and ylist are randomly sorted per subsquare, like in latinhypercube
    for (i=0; i<N_sqrt; i++)
    {
	for(j=0; j<N_sqrt; j++) 
	{ 
	    number = (((long) genrand_int31() % (N_sqrt-j)) + j), old = xlist[i][j];
	    xlist[i][j] = xlist[i][number];
	    xlist[i][number] =  old;
	}
    }
    for (i=0; i<N_sqrt; i++)
    {
	for(j=0; j<N_sqrt; j++) 
	{ 
	    number = (((long) genrand_int31() % (N_sqrt-j)) + j), old = ylist[i][j];
	    ylist[i][j] = ylist[i][number];
	    ylist[i][number] =  old;
	}
    }
    
    //Loop over subsquares
    for (i=0; i<N_sqrt; i++)
    {
	//Loop over columns/rows within a subsquare
	for (j=0; j<N_sqrt; j++)
	{
	    //Generate random number in subsquare (i,j) with column j and row i
	    x = (xlist[i][j] + genrand_real2()) * x_Scale + MinRe;
	    y = (ylist[j][i] + genrand_real2()) * y_Scale + MinIm;
	    Z_re = x, Z_im = y;

	    //Loop over iterations
	    for(n=0; n<=MaxIterations; n++)
	    {
		Z_re2 = Z_re*Z_re, Z_im2 = Z_im*Z_im;
		if((Z_re2 + Z_im2) > 4) { break; }
		Area_OR_It[n] += 1; //count samples that are still in Mandelbrot set at iteration i
		Z_im = 2*Z_re*Z_im + y;
		Z_re = Z_re2 - Z_im2 + x;
	    }
	}
    }
    //Calculate area as function of iteration, add to mean area (squared)
    for(i=0; i<(MaxIterations+1); i++) 
    { 
	Area_OR_It[i] *= (Surface/N); 
	Area_OR_Mean[i] += Area_OR_It[i];
	Area_OR_Sq_Mean[i] += Area_OR_It[i]*Area_OR_It[i];
    }
}

//Main body
int main()
{
    cout.precision(6);
    
    //File writing
    //ofstream file;
    //file.open("N_10E6_5.dat");
     
    //Initialize arrays for the mean area, squared mean area and the error
    vector<float> Area_RS_Mean (MaxIterations+1), Area_RS_Sq_Mean (MaxIterations+1), Area_RS_Error (MaxIterations+1);
    vector<float> Area_LH_Mean (MaxIterations+1), Area_LH_Sq_Mean (MaxIterations+1), Area_LH_Error (MaxIterations+1);
    vector<float> Area_OR_Mean (MaxIterations+1), Area_OR_Sq_Mean (MaxIterations+1), Area_OR_Error (MaxIterations+1);
   
    //Loop over runs
    for(int m=0; m<N_runs; m++)
    {
	if(m%(N_runs/100)==0) cout << (float) m/N_runs * 100 << "% done" << endl;

	//Set seed 
	init_genrand(m+time(NULL));
		
	//Initialize arrays which give area per run, index = iteration
	vector<float> Area_RS_It (MaxIterations+1),Area_LH_It (MaxIterations+1),Area_OR_It (MaxIterations+1);
	
	//Do the sampling and fill the arrays above
	if(RS==true) random(Area_RS_It,Area_RS_Mean,Area_RS_Sq_Mean);
	if(LH==true) latinhypercube(Area_LH_It,Area_LH_Mean,Area_LH_Sq_Mean);
	if(OR==true) orthogonal(Area_OR_It,Area_OR_Mean,Area_OR_Sq_Mean);
    }

    //After doing all runs again loop over all iterations
    if(RS==true)
    {
	for(int i=0; i<MaxIterations+1; i++)
	{       
	    Area_RS_Mean[i] /= N_runs; Area_RS_Sq_Mean[i] /= N_runs;
	    Area_RS_Error[i] = (2.58/sqrt(N_runs)) * (sqrt(N_runs/(N_runs-1))) * sqrt(Area_RS_Sq_Mean[i]-(Area_RS_Mean[i]*Area_RS_Mean[i]));
	    //file << i << " " << Area_RS_Mean[i] << " " <<  Area_RS_Error[i] << "\n";
	    if(i%20==0) { cout << "RS Iteration " << i << " Area = " << Area_RS_Mean[i] << " Error = " << Area_RS_Error[i] << endl; }
	}
    }
    if(LH==true)
    {
	for(int i=0; i<MaxIterations+1; i++)
	{
	    Area_LH_Mean[i] /= N_runs; Area_LH_Sq_Mean[i] /= N_runs;
	    Area_LH_Error[i] = (2.58/sqrt(N_runs)) * (sqrt(N_runs/(N_runs-1))) * sqrt(Area_LH_Sq_Mean[i]-(Area_LH_Mean[i]*Area_LH_Mean[i]));
	    //file << i << " " << Area_LH_Mean[i] << " " <<  Area_LH_Error[i] << "\n";
	    if(i%20==0) { cout << "LH Iteration " << i << " Area = " << Area_LH_Mean[i] << " Error = " << Area_LH_Error[i] << endl; }
	}
    }
    if(OR==true)
    {
	for(int i=0; i<MaxIterations+1; i++)
	{
	    Area_OR_Mean[i] /= N_runs; Area_OR_Sq_Mean[i] /= N_runs;
	    Area_OR_Error[i] = (2.58/sqrt(N_runs)) * (sqrt(N_runs/(N_runs-1))) * sqrt(Area_OR_Sq_Mean[i]-(Area_OR_Mean[i]*Area_OR_Mean[i]));
	    //file << i << " " << Area_OR_Mean[i] << " " <<  Area_OR_Error[i] << "\n";
	    if(i%20==0) { cout << "OR Iteration " << i << " Area = " << Area_OR_Mean[i] << " Error = " << Area_OR_Error[i]<< endl; }
	}
    } 
   //file.close();
}
