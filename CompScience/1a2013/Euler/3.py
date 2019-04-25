number = 600851475143

x = 2
while number > 1:
    if number%x == 0:
        number = number/x
        print x
    x = x + 1
        

print 'The largest prime factor of 600851475143 is',x-1
