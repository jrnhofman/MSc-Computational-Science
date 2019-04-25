import numpy as np
import timeit
import csv
from scipy.linalg import lu
time=0
n=100
f = open('/home/jhofman/Desktop/3_3.txt','w')
while n < 1000:
    for i in range(0,50):
        t = timeit.Timer("(P,L,U) = lu(A)","import numpy as np; from scipy.linalg import lu; A = np.matrix(np.random.randn("+str(n)+","+str(n)+"))")
        time += t.timeit(1)
        if i == 49:
            time = time/50
            print time
            s = str(n)+' '+str(time)+' '
            f.write(s)
            
    n+=100

