hundred = 7

#1 through 19
first = [3,3,5,4,4,3,5,5,4,3,6,6,8,8,7,7,9,8,8]

# 10 through 90
second = [3,6,6,5,5,5,7,6,6]

#1 through 19
number = sum(first)

#20 through 99
number += 8*sum(first[:9]) + 10*(sum(second)-second[0])

#100 through 999
number += 9*sum(first) + 9*(8*sum(first[:9]) + 10*(sum(second)-second[0])) + 900*hundred + 891*3 + 100*sum(first[:9])

#1000
number += 11

print "The sum of all the letters from 1 till 1000 = ",number
