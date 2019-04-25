#include <stdio.h>
#include <math.h>

int i;
double h,j;
double deriv[50];

main(){
	FILE *f;
	f = fopen("deriv2.txt","w");
	j=1;
	for(i=1;i<50;i++){
		j*=2;
		h=1/(double)j;
		deriv[i]=-((sin(1.+h)-sin(1.-h))/(2*h))+cos(1);
		if (deriv[i]<0) deriv[i]*=-1;
		fprintf(f,"%e %.16f\n",h,deriv[i]);
		printf("%e %.16f\n",h,deriv[i]);
	}
	fclose(f);
}
		
