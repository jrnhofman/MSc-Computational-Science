import numpy as np
import scipy.linalg
from scipy.linalg import solve


A = np.matrix([[1.,0.,0.,0.,0.,0.],[1.,0.5,0.25,0.125,0.06125,0.030625],[1.,1.,1.,1.,1.,1.],[1.,6.,36.,216.,1296.,7776.],[1.,7.,49.,343.,2401.,16807.],[1.,9.,81.,729.,6561.,59049.]])
y = np.array([0.0,1.6,2.0,2.0,1.5,0.0])

x = solve(A,y)
print x
