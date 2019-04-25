import random
file = open('/home/jhofman/Desktop/Project/basicfb.txt','w')

S0,E0 = 1.,1.
kp,km,r,rm = random.random(),random.random(),random.random(),random.random()
KM = (km+r)/kp
d = 0.01

S = S0
E = E0
C = 0.
P = 0.

def f(S,E,C,P):
    y = -kp*E*S + km*C + rm*E*P
    return y
def g(S,E,C,P):
    y = -kp*E*S + (km + r)*C
    return y
def h(S,E,C,P):
    y = kp*E*S - (km + r)*C
    return y
def i(S,E,C,P):
    y = r*C - rm*P*E
    return y

t = 0
tmax = 10
while t < tmax:
    file.write(str(t)+' '+str(S)+' '+str(E)+' '+str(C)+' '+str(P)+'\n')
    print "t = ",t," S = ",S," E = ",E," C = ",C," P = ",P
    k1 = d*f(S,E,C,P)
    l1 = d*g(S,E,C,P)
    m1 = d*h(S,E,C,P)
    n1 = d*i(S,E,C,P)
    #print k1,l1,m1,n1
    k2 = d*f(S+k1/2,E+l1/2,C+m1/2,P+n1/2)
    l2 = d*g(S+k1/2,E+l1/2,C+m1/2,P+n1/2)
    m2 = d*h(S+k1/2,E+l1/2,C+m1/2,P+n1/2)
    n2 = d*i(S+k1/2,E+l1/2,C+m1/2,P+n1/2)
    #print k2,l2,m2,n2
    k3 = d*f(S+k2/2,E+l2/2,C+m2/2,P+n2/2)
    l3 = d*g(S+k2/2,E+l2/2,C+m2/2,P+n2/2)
    m3 = d*h(S+k2/2,E+l2/2,C+m2/2,P+n2/2)
    n3 = d*i(S+k2/2,E+l2/2,C+m2/2,P+n2/2)
    #print k3,l3,m3,n3
    k4 = d*f(S+k3,E+l3,C+m3,P+n3)
    l4 = d*g(S+k3,E+l3,C+m3,P+n3)
    m4 = d*h(S+k3,E+l3,C+m3,P+n3)
    n4 = d*i(S+k3,E+l3,C+m3,P+n3)
    #print k4,l4,m4,n4
    k = (1./6.)*(k1+2*k2+2*k3+k4)
    l = (1./6.)*(l1+2*l2+2*l3+l4)
    m = (1./6.)*(m1+2*m2+2*m3+m4)
    n = (1./6.)*(n1+2*n2+2*n3+n4)
    #print k,l,m,n
    S += k
    E += l
    C += m
    P += n
    t += d

file.close()
