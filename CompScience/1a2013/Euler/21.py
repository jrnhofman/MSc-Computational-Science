from math import sqrt

def divisors(i):
    sum = 1
    x = 2
    while x < sqrt(i):
        if i%x==0:
            sum += x
            sum += i/x       
        x += 1
    if i%sqrt(i)==0:
        sum += int(sqrt(i))
    return sum

xmax = 10000
amicable = 0

divisorsum = [0 for x in range(xmax)]
for i in range(xmax-1,0,-1):
    divisorsum[i] = divisors(i+1)

for i in range(0,xmax):
    for j in range(i+1,xmax):
        a = divisorsum[i]
        b = divisorsum[j]
        if a==j+1 and b==i+1 and i+1!=j+1:
            print a,b,i,j
            amicable += i+1
            amicable += j+1
print 'The sum of amicable numbers up to 10000 is',amicable
