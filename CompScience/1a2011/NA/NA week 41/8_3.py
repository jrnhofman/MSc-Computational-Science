from math import cos
from math import sin
from math import atan

def midpoint(xlist,f):
    sum = 0
    for i in range(1,len(xlist)):
        sum += f((xlist[i-1]+xlist[i])/2)
    sum *= 2/float(n)
    return sum

def trapezoid(xlist,f):
    sum = 0
    for i in range(1,len(xlist)):
        sum += f(xlist[i-1]) + f(xlist[i])
    sum *= 1/float(n)
    return sum

def simpson(xlist,f):
    sum = 0
    for i in range(1,len(xlist)):
        sum += f(xlist[i-1]) + 4*f((xlist[i-1]+xlist[i])/2) + f(xlist[i])
    sum *= 1/(3.*float(n))
    return sum

def fa(x):
    y = cos(x)
    return y
def fb(x):
    y = 1/(1+100*x*x)
    return y

sola = 2*sin(1)
solb = 0.2*atan(10) 
nlist = [n*2 for n in range(1,11)]
for n in nlist:
    xlist = [-1 + 2*x/float(n) for x in range(0,n+1)] #points (= #intervals+1)
    ma = midpoint(xlist,fa)
    mb = midpoint(xlist,fb)
    ma_error = abs(ma-sola)
    mb_error = abs(mb-solb)
    ta = trapezoid(xlist,fa)
    tb = trapezoid(xlist,fb)
    ta_error = abs(ta-sola)
    tb_error = abs(tb-solb)
    sa = simpson(xlist,fa)
    sb = simpson(xlist,fb)
    sa_error = abs(sa-sola)
    sb_error = abs(sb-solb)
    #print "n = ",n," M(f) = ",ma," T(f) = ",ta," S(f) = ",sa
    print "n = ",n," error M(f) = ",mb_error," error T(f) = ",tb_error," error S(f) =",sb_error
    if mb_error == 0:
        break
