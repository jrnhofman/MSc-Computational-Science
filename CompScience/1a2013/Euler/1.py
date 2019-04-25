sum = 0
for x in range(1,1000):
    if x%3 == 0:
        sum += x
    elif x%5==0:
        sum += x

print 'The sum of all numbers divisible by 3 and 5 below 1000 is',sum
