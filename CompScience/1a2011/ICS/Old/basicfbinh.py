import random
file = open('/home/jhofman/Desktop/Project/basicfbinh.txt','w')

S0,E0,I0,J0 = 1.,1.,1.,1.
kp,km,r,rm,qp,qm = random.random(),random.random(),random.random(),random.random(),random.random(),random.random()
KM = (km+r)/kp
d = 0.01

S = S0
E = E0
I = I0
C = 0.
J = J0
P = 0.

def f(S,E,C,P,I,J):
    y = -kp*E*I*S + km*C*I
    return y
def g(S,E,C,P,I,J):
    y = qp*C*S - qm*E*I*S + r*C*I - rm*E*P*I + km*C*I - kp*E*I*S
    return y
def h(S,E,C,P,I,J):
    y = kp*E*I*S + rm*E*P*I - r*C*I - km*C*I
    return y
def i(S,E,C,P,I,J):
    y = r*C - rm*E*P*I
    return y
def o(S,E,C,P,I,J):
    y = qp*J*S - qm*E*I*S
    return y
def p(S,E,C,P,I,J):
    y = qm*E*I*S - qp*J*S
    return y

t = 0
tmax = 25
while t < tmax:
    file.write(str(t)+' '+str(S)+' '+str(E)+' '+str(C)+' '+str(P)+' '+str(I)+' '+str(J)+'\n')
    print "t = ",t," S = ",S," E = ",E," C = ",C," P = ",P," I = ",I," J = ",J
    k1 = d*f(S,E,C,P,I,J)
    l1 = d*g(S,E,C,P,I,J)
    m1 = d*h(S,E,C,P,I,J)
    n1 = d*i(S,E,C,P,I,J)
    u1 = d*o(S,E,C,P,I,J)
    s1 = d*p(S,E,C,P,I,J)
    #print k1,l1,m1,n1,u1,s1
    k2 = d*f(S+k1/2,E+l1/2,C+m1/2,P+n1/2,I+u1/2,J+s1/2)
    l2 = d*g(S+k1/2,E+l1/2,C+m1/2,P+n1/2,I+u1/2,J+s1/2)
    m2 = d*h(S+k1/2,E+l1/2,C+m1/2,P+n1/2,I+u1/2,J+s1/2)
    n2 = d*i(S+k1/2,E+l1/2,C+m1/2,P+n1/2,I+u1/2,J+s1/2)
    u2 = d*o(S+k1/2,E+l1/2,C+m1/2,P+n1/2,I+u1/2,J+s1/2)
    s2 = d*p(S+k1/2,E+l1/2,C+m1/2,P+n1/2,I+u1/2,J+s1/2)
    #print k2,l2,m2,n2,u1,s1
    k3 = d*f(S+k2/2,E+l2/2,C+m2/2,P+n2/2,I+u2/2,J+s2/2)
    l3 = d*g(S+k2/2,E+l2/2,C+m2/2,P+n2/2,I+u2/2,J+s2/2)
    m3 = d*h(S+k2/2,E+l2/2,C+m2/2,P+n2/2,I+u2/2,J+s2/2)
    n3 = d*i(S+k2/2,E+l2/2,C+m2/2,P+n2/2,I+u2/2,J+s2/2)
    u3 = d*o(S+k2/2,E+l2/2,C+m2/2,P+n2/2,I+u2/2,J+s2/2)
    s3 = d*p(S+k2/2,E+l2/2,C+m2/2,P+n2/2,I+u2/2,J+s2/2)
    #print k3,l3,m3,n3
    k4 = d*f(S+k3,E+l3,C+m3,P+n3,I+u3,J+s3)
    l4 = d*g(S+k3,E+l3,C+m3,P+n3,I+u3,J+s3)
    m4 = d*h(S+k3,E+l3,C+m3,P+n3,I+u3,J+s3)
    n4 = d*i(S+k3,E+l3,C+m3,P+n3,I+u3,J+s3)
    u4 = d*o(S+k3,E+l3,C+m3,P+n3,I+u3,J+s3)
    s4 = d*p(S+k3,E+l3,C+m3,P+n3,I+u3,J+s3)
    #print k4,l4,m4,n4
    k = (1./6.)*(k1+2*k2+2*k3+k4)
    l = (1./6.)*(l1+2*l2+2*l3+l4)
    m = (1./6.)*(m1+2*m2+2*m3+m4)
    n = (1./6.)*(n1+2*n2+2*n3+n4)
    u = (1./6.)*(u1+2*u2+2*u3+u4)
    s = (1./6.)*(s1+2*s2+2*s3+s4)
    #print k,l,m,n
    S += k
    E += l
    C += m
    P += n
    I += u
    J += s
    t += d

file.close()
