from math import floor
import time

start = time.time()

primelist = [2]
n = 3
while n < 200000:
    check = True
    for x in primelist:
        if x > floor(n**0.5):
            break
        if n%x == 0:
            n = n + 2
            check = False
            break
    if check:
        primelist.append(n)
        n = n + 2

end = time.time()
print 'Time elapsed:',end-start

print 'The sum of all primes beneath 2000000 is',sum(primelist)


