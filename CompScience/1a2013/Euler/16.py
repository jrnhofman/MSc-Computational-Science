print str(2**1000)[0]

sum = 0
print len(str(2**1000))
for x in range(len(str(2**1000))):
    sum += int(str(2**1000)[x])
print 'The sum of the digits of 2**1000 is',sum
