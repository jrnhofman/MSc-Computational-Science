import numpy as np
import scipy.linalg
from scipy.linalg import qr

A42 = np.matrix([[2,3,2],[10,3,4],[3,6,1]])
A43 = np.matrix([[6,2,1],[2,3,1],[1,1,1]])

n=10
for j in range(0,n):
    sigma = float(A43[2,2])
    (Q,R) = qr(A43-sigma*np.identity(3))
    A43 = np.dot(R,Q) + sigma*np.identity(3)
    sigma2 = float(A43[2,2])
    print sigma
    print sigma2/sigma
