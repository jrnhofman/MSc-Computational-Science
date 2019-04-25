import numpy as np
import csv
from numpy.linalg import cholesky
t = np.array([0.0,1.0,2.0,3.0,4.0,5.0])

n=1;
f = csv.writer(open('/home/jhofman/Desktop/3_1.txt', 'w'), delimiter=' ')
while n < 7:
    k = 0
    b = np.array([[1.0],[2.7],[5.8],[6.6],[7.5],[9.9]])
    A = np.matrix(np.zeros([6,n]))
    y = np.array(np.zeros([n]))
    x = np.array(np.zeros([n]))
    for i in range(0,6):
        for j in range(0,n):
            A[i,j]=t[i]**j
    b = np.dot(np.transpose(A),b)
    A = np.dot(np.transpose(A),A)
    L = cholesky(A)

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
        
    string=(n,x)
    s = str(string)
    f.writerow(x)
    n+=1

