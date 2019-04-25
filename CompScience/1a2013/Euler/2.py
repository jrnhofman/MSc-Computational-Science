from time import time,clock

start = time()

a = 1
b = 2
c = 0

sum = b
while c < 4000000:
    c = a + b
    if c&1 == 0:
        sum += c
    a = b
    b = c

end = (time() - start)
print 'Elapsed time',end

print 'The sum of all even Fibonacci numbers under 4 million is',sum
