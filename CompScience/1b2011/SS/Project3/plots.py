import matplotlib.pyplot as plt
import matplotlib.mlab as lab
from numpy import transpose

X1,X2,X3 = transpose(lab.load("ElevatorSample.dat",unpack=True))
Y1,Y2,Y3 = transpose(lab.load("ElevatorSamplePass.dat",unpack=True))
Z1,Z2,Z3 = transpose(lab.load("ElevatorSamplePassMult.dat",unpack=True))
plt.figure(1)
plt.errorbar(X1,X2,X3,fmt='b.',label='case 1')
plt.errorbar(Y1,Y2,Y3,fmt='r.',label='case 2')
plt.errorbar(Z1,Z2,Z3,fmt='y.',label='case 3')
plt.xlabel("Mean time between arrivals")
plt.ylabel("Average sample difference")
plt.legend(loc=4)


X1,X2,X3 = transpose(lab.load("Elevator.dat",unpack=True))
Y1,Y2,Y3 = transpose(lab.load("ElevatorPass.dat",unpack=True))
Z1,Z2,Z3 = transpose(lab.load("ElevatorPassMult.dat",unpack=True))
X1 = list(X1)
X2 = list(X2)
X3 = list(X3)
Y1 = list(Y1)
Y2 = list(Y2)
Y3 = list(Y3)
Z1 = list(Z1)
Z2 = list(Z2)
Z3 = list(Z3)
del X1[0:15]
del X2[0:15]
del X3[0:15]
del Y1[0:15]
del Y2[0:15]
del Y3[0:15]
del Z1[0:15]
del Z2[0:15]
del Z3[0:15]
plt.figure(2)
plt.errorbar(X1,X2,X3,fmt='b.',label='case 1')
plt.errorbar(Y1,Y2,Y3,fmt='r.',label='case 2')
plt.errorbar(Z1,Z2,Z3,fmt='y.',label='case 3')
plt.xlabel("Mean time between arrivals")
plt.ylabel("Average journey time")
plt.legend(loc=1)
plt.show()
