#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <stdlib.h>
using namespace std;

//Resolution of grid
int ImageHeight = 1000;
int ImageWidth = 1000;
int MaxIterations = 100;

//Parameters for dimension of image
float MinRe = -2.0;
float MaxRe = 1.0;
float MinIm = -1.2;
float MaxIm = MinIm+(MaxRe-MinRe)*ImageHeight/ImageWidth;
float Re_factor = (MaxRe-MinRe)/(ImageWidth-1);
float Im_factor = (MaxIm-MinIm)/(ImageHeight-1);

int main()
{
    //File writing
    ofstream file;
    file.open ("graphics.dat"); 

    //Loop over all pixels
    for(int y=0; y<ImageHeight; y++)
    {
	float c_im = MaxIm - y*Im_factor;
	if(y!=0){ file << "\n"; }
	
	for(int x=0; x<ImageWidth; x++)
	{
	    float c_re = MinRe + x*Re_factor;
	    float Z_re = c_re, Z_im = c_im;
	    
	    //check if (x,y) is in Mandelbrot set and save iteration number when it jumps out
	    //this routine saves 0 when a pixel is always in the set (regardless of iteration)
	    //and also saves 0 if it is not in the set after MaxIterations, otherwise it returns
	    //the iteration
	    for(int n=0; n<(MaxIterations+1); n++)
	    {
		float Z_re2 = Z_re*Z_re, Z_im2 = Z_im*Z_im;
		if(Z_re2 + Z_im2 > 4)
		{
		    file << float(n)/MaxIterations << " ";
		    break;
		}
		if(n == MaxIterations && (Z_re2 + Z_im2) <= 4) { file << 0 << " "; }
		Z_im = 2*Z_re*Z_im + c_im;
		Z_re = Z_re2 - Z_im2 + c_re;
	    }
	}
    }
    file.close();
}
