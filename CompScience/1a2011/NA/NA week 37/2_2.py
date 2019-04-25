import numpy as np
from scipy.linalg import lu

A = np. matrix('2. 4. -2.;4. 9. -3.;-2. -1. 7.')
b = np.array([4,8,-6])
c = np.array([4,8,-6])
u = np.array([1,0,0])
v = np.array([0,2,0])
x = np.array(np.zeros([3]))
y = np.array(np.zeros([3]))
z = np.array(np.zeros([3]))
q = np.array(np.zeros([3]))
n = 3

#PLU factorization of A
(P,L,U) = lu(A)

#Solve PLy = b
L = np.dot(P,L)
y[0] = b[1]/L[1,0]
y[1] = (b[2] - L[2,0]*y[0])/L[2,1]
y[2] = (b[0] - L[0,0]*y[0] - L[0,1]*y[1])/L[0,2]

#Solve Ux = y
x[2] = y[2]/U[2,2]
x[1] = (y[1] - U[1,2]*x[2])/U[1,1]
x[0] = (y[0] - U[0,1]*x[1] - U[0,2]*x[2])/U[0,0]

print x
#2.2c Take u = (1,0,0)^T and v = (0,2,0), then solve Az = PLUz = u
#Solve PLq = u
q[0] = u[1]/L[1,0]
q[1] = (u[2] - L[2,0]*q[0])/L[2,1]
q[2] = (u[0] - L[0,0]*q[0] - L[0,1]*q[1])/L[0,2]

#Solve Uz = q
z[2] = q[2]/U[2,2]
z[1] = (q[1] - U[1,2]*z[2])/U[1,1]
z[0] = (q[0] - U[0,1]*z[1] - U[0,2]*z[2])/U[0,0]

x = x + np.dot(v,x)/(1-np.dot(v,z)) * z

print x




# while k < n:
#    if A[k,k] is not 0:
#        for i in range(k+1,n+1):
#            m[i,k] = A[i,k]/A[k,k]
#    for i in range (k+1,n+1):
#        b[i] = b[i] - m[i,k]*b[k]
#        for j in range (k,n+1):
#            A[i,j] = A[i,j] - m[i,k]*A[k,j]
#    k +=1
#print A
#print b
