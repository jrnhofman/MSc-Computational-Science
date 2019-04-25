import string
a = "".join(map(string.strip,open('7.dat').readlines()))
zipped = zip(a[1:],a[2:],a[3:],a[4:],a[5:])
mapped = map(lambda (x1,x2,x3,x4,x5): int(x1)*int(x2)*int(x3)*int(x4)*int(x5),zipped)
print 'The maximum of the product of five consecutive numbers is',max(mapped)
