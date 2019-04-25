# We have a**2 + b**2 = c**2 and a + b + c = 1000, 
# substituting the second into the first gives after 
# solving that b = (2000a-10**6)/(2a-2000), which excludes 
# a > 499 since b has to be positive. b has to be integer 
# hence the quotient should be an integer 
# (so use modulo operation) 

for a in range(1,500):
    if (2000*a-10**6)%(2*a-2000)==0:
        b = (2000*a-10**6)/(2*a-2000)
        c = 1000 - a - b
        break
print 'The Pythagorean triplet for which a + b + c = 1000 is',a,b,c
print 'Their product is',a*b*c
