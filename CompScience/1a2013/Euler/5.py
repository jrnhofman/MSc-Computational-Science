# Crappy solution, takes ages

# number = 2520

# i = 2
# while i < 21:
#    if number%i == 0:
#        i = i + 1
#    else:
#        number += 20
#        print number
#        i = 2

# print 'The smallest number for which 1..20 is a divisor is',number

def lcm(a,b):
    gcd,tmp = a,b
    while tmp != 0:
        gcd,tmp = tmp,gcd%tmp
    return a*b/gcd

numbers = [x for x in range(1,21)]
number = reduce(lcm,numbers)

print 'The smallest number for which 1..20 is a divisor is',number
