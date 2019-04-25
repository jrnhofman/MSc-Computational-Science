#include <signal.h>
#include <unistd.h>
#include <stdio.h>

#define TIME 15

static volatile int tag = 1;

void handler(int unused) {
 tag = 0;
}

int main(void) {
 long long i = 0;
 double sign = 1.0;
 double odd  = 1.0;
 double sum  = 0.0;
 double res  = 0.0;

 signal(SIGALRM, handler);
 alarm(TIME);

 while (tag) {
   double frac = sign / odd;
   sum += frac;
   res = 4.0 * sum;
   sign *= -1.0;
   odd += 2.0;
   i += 1;
 }

 printf("%lld %.15lf %lld\n", 
        i, res, (i*5)/TIME);
 return 0;
}
