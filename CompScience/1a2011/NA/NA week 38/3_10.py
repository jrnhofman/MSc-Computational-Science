import numpy as np
import random
import scipy.linalg
from scipy.linalg import inv

#Generate random matrix A
m = random.randint(3,5)
n = random.randint(m,5)
A = np.random.randn(n,m)
print A

#QR decomposition of A
(Q,R) = np.linalg.qr(A)
R_t = np.transpose(R)
print R

#Compute inverse of R
R_inv = np.matrix(np.zeros([m,m]))
for i in range(0,m):
    R_inv[i,i] = 1/R[i,i]
for i in range(0,m):
    for j in range(i+1,m):
        R_inv[i,j] = 0
        for k in range(i,j):
            R_inv[i,j] -= R_inv[i,k]*R[k,j]
        R_inv[i,j] /= R[j,j]
R_t_inv = np.transpose(R_inv)

#Compute covariance matrix and compare with (A^T A)^-1
print np.dot(R_inv,R_t_inv)
print inv(np.dot(np.transpose(A),A))
    

