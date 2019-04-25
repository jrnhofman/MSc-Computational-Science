
import numpy as np
import scipy.linalg
from scipy.linalg import lu,eig

A = np.matrix([[6,2,1],[2,3,1],[1,1,1]])
x = np.transpose([1,1,1])

#function to compute infinity-norm
def norm(z):
    sum = 0
    for j in range(0,len(z)):
        sum += abs(float(z[[j]]))
    return sum


#function to solve LU x = b, see computer exercise 2.6
def solve(P,L,U,b):
    q = np.array(np.zeros([3]))  
    z = np.array(np.zeros([3]))
    b = np.dot(np.transpose(P),b)
    k = 0
    while k < 3:
        q[k] = b[k]
        for l in range(0,k):
            q[k] -= L[k,l]*q[l]
        q[k] = q[k]/L[k,k]
        k += 1

    k = 2
    while k > -1:
        z[k] = q[k]
        for l in range(k+1,3):
            z[k] -= U[k,l]*z[l]
        z[k] = z[k]/U[k,k]
        k -= 1
    return z


#PLU decomposition of shifted A
(P,L,U) = lu(A-2*np.identity(3))

#inverse iteration
n = 10
for j in range(0,n):
    y = solve(P,L,U,x)
    x = y/norm(y)
    print x
    print 1/norm(y)

#use library routine to calculate eigenvalues and -vectors
(w,v) = eig(A)
print w
print v

#Rayleigh quotient iteration
x = np.transpose([1,1,1])
m = 10

for j in range(0,m):
    sigma = np.dot(np.dot(np.transpose(x),A),x)/np.dot(np.transpose(x),x)
    (P,L,U) = lu(A-float(sigma)*np.identity(3))
    y = solve(P,L,U,x)
    x = y/norm(y)
    print x
    print sigma
