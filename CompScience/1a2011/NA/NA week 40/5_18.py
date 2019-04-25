import numpy as np
from math import sin
from math import cos
from math import exp
from scipy.linalg import solve

def g(x,y):
    gx = (x+3)*(y*y-7)+18
    gy = sin(y*exp(x)-1)
    return np.array([gx,gy])

def Jac(x,y):
    j1 = y**3 - 7
    j2 = 3*(x+3)*y*y
    j3 = exp(x)*y*cos(exp(x)*y-1)
    j4 = exp(x)*cos(exp(x)*y-1)
    return np.matrix([[j1,j2],[j3,j4]])

x0,y0 = -0.5,1.4
n = 20
for k in range(0,n):
    print x0,y0
    print "k = ",k," error = ",x0+y0-1
    s = solve(Jac(x0,y0),-g(x0,y0))
    x0 += s[0]
    y0 += s[1]

x0,y0 = -0.5,1.4
Jaco = Jac(x0,y0) 
n = 100
for k in range(0,n):
    print x0,y0
    print "k = ",k," error = ",x0+y0-1
    s = np.transpose(solve(Jaco,-g(x0,y0)))
    x0 += s[0]
    y0 += s[1]
    z = g(x0,y0)-g(x0-s[0],y0-s[1])
    Jaco = Jaco + (np.outer((z - np.dot(Jaco,s)),s))/(np.dot(s,s))
     


