from math import sqrt
from math import log

def g1(x):
    y = (x*x+2.)/3.
    return y
def g2(x):
    y = sqrt(3.*x-2.)
    return y
def g3(x):
    y = 3.-2./x
    return y
def g4(x):
    y = (x*x-2.)/(x*2.-3.)
    return y

x1,x2,x3,x4 = 3.,3.,3.,3.
c2,c3 = 0.75,0.5
n = 20
x = 2
for j in range(0,n):
    err2k = abs(x-x2)
    err3k = abs(x-x3)
    x1 = g1(x1)
    x2 = g2(x2)
    x3 = g3(x3)
    x4 = g4(x4)
    err2kp = abs(x-x2)
    err3kp = abs(x-x3)
    if j==n-1:
        r2 = log(err2kp/c2,err2k)
        r3 = log(err3kp/c3,err3k)
        print "r2 = ",r2," r3 = ",r3
        print "x1 = ",x1," x2 = ",x2," x3 = ",x3," x4 = ",x4
 

    
