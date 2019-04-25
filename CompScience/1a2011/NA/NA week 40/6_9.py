import numpy as np
from scipy.linalg import solve

f = open('/home/jhofman/Desktop/NA/NA week 40/path3.txt','w')

def gradg(x,y):
    gx = 200*(y-x*x)*-2*x -2*(1-x)
    gy = 200*(y-x*x)
    return np.array([gx,gy])

def Hes(x,y):
    j1 = 1200*x*x-400*y+2
    j2 = -400*x
    j3 = -400*x
    j4 = 200
    return np.matrix([[j1,j2],[j3,j4]])

x0,y0 = 2.0,1.0
n = 10
for k in range(0,n):
    print x0,y0
    f.write(str(x0)+' '+str(y0)+'\n')
    s = solve(Hes(x0,y0),-gradg(x0,y0))
    x0 += s[0]
    y0 += s[1]
    
