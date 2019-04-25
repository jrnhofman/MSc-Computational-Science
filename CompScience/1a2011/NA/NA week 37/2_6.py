import numpy as np
from numpy.linalg import cholesky

f = open('/home/jhofman/Desktop/2_6.txt','w')

n=3
while n<14:
    k=0
    A = np.array(np.zeros([n,n]))
    y = np.array(np.zeros([n]))
    x = np.array(np.zeros([n]))
    #generate Hilbert matrix of order n
    for i in range (0,n):
        for j in range (0,n):
            A[i,j] = 1/(-1+(float(i)+1)+(float(j)+1))

    xsol = np.array(np.ones([n]))
    b = np.dot(A,xsol)
    L  = cholesky(A)

#solve Ly = b
    while k < n:
        y[k] = b[k]
        for l in range(0,k):
            y[k] -= L[k,l]*y[l]
        y[k] = y[k]/L[k,k]
        k += 1

#solve L^Tx = y
    U = np.transpose(L)
    k = n -1
    while k > -1:
        x[k] = y[k]
        for l in range(k+1,n):
            x[k] -= U[k,l]*x[l]
        x[k] = x[k]/U[k,k]
        k -= 1
    
#compute errors
    delta = x - xsol
    r = b - np.dot(A,x)
    residual,error = 0,0
    for j in range (0,n):
        residual+=abs(r[j])
        error+=abs(delta[j])
    f.write(str(n)+' '+str(error)+' ')
    n +=1
f.closed
