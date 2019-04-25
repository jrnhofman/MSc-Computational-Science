sumsq,sqsum = 0,0

for x in range(1,101):
    sumsq += x*x
    sqsum += x

print 'The squared sum minus the sum of squares for 1..100 equals',sqsum*sqsum-sumsq
