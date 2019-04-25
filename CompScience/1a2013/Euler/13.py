a = open('13.dat').readlines()

int_sum = 0
for x in range(len(a)):
    int_sum += int(a[x])
int_sum = str(int_sum)[:10]
print 'The first 10 digits of the sum are',int_sum
