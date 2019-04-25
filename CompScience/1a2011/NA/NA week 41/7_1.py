import random
import numpy as np

def eval(n,x,t,deriv=False,inte=False,a=0,b=0):
    if deriv==False and inte==False:
        value  = x[-1]
        for j in range(2,n+2):
            value *= t
            value += x[-j]
        return value
    if deriv == True:
        y = [0]*(len(x)-1)
        n = n-1
        for i in range(0,len(x)-1):
            y[i] = x[i+1]*(i+1)
        value = y[-1]
        for j in range(2,n+2):
            value *= t
            value += y[-j]
        return value
    if inte==True and a!=0 and b!=0:
        y = [0]*(len(x)+1)
        n = n + 1
        for i in range(1,len(x)+1):
            y[i] = float(x[i-1])/float(i)
        value1,value2 = y[-1],y[-1]
        for j in range(2,n+2):
            value1 *= a
            value2 *= b
            value1 += y[-j]
            value2 += y[-j]
        return value2-value1
    else: 
        return 'wrong input!'

n = random.randint(0,10)
x = np.random.randint(0,10,n+1)
t = random.randint(0,10)
value = eval(2,[2,-6,2],2,False,True,2,4)
print value

