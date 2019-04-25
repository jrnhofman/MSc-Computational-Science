import numpy as np
import scipy.linalg
from scipy.linalg import inv
from numpy import outer
from scipy.linalg import eig

#Define matrix and the starting vectors
A = np.matrix([[2,3,2],[10,3,4],[3,6,1]])
x = np.array([[0],[0],[1]])
r = np.array([[0],[0],[1]]) 

#function to compute infinity-norm
def norm(z):
    sum = 0
    for j in range(0,len(z)):
        sum += float(z[[j]])
    return sum

n=1000

#(normed) power iteration
for j in range(0,n):
    y = np.dot(A,x)
    x = y/norm(y)
    if j == n-1:
        print x #eigenvector
        print norm(y) #eigenvalue

#make new matrix with largest eigenvalue removed and apply power iteration to that matrix
B = A - np.outer(np.transpose(x),np.array([10,10,13]))

for j in range(0,n):
    y = np.dot(B,r)
    r = y/norm(y)
    if j == n-1:
        print r #eigenvector
        print norm(y) #eigenvalue

#use library routine to calculate eigenvalues and -vectors
(w,v) = eig(A)
print w
print v
