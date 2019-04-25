a = open('18.dat').readlines()
for i in range(len(a)):
    a[i] = [int(x) for x in a[i].split()]

for i in range(len(a)-2,-1,-1):
    for j in range(i+1):
        a[i][j] += max(a[i+1][j],a[i+1][j+1])
print 'The maximum sum in the triangle is',a[0][0]
