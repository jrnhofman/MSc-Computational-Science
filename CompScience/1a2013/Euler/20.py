from math import factorial

sum = 0
for i in range(len(str(factorial(100)))):
    sum += int(str(factorial(100))[i])
print 'The sum of the digits of 100! is',sum

