number = []
for x in range(10**6):
    number.append(0)

for x in range(10**6-1,0,-1):
    if number[x]==0:
        add = 1
        y = x
        while True:
            if y==1:
                number[x] = add
                break
            elif y%2 == 0:
                y /= 2
                if y < 10**6:
                    number[y] = -1
                add += 1
            else:
                y = 3*y + 1
                if y < 10**6:
                    number[y] = -1
                add += 1
print 'The maximum chain length belongs to the number',number.index(max(number))

# Faster solution from internet, this one uses storing the whole chain efficiently (because of the returns of the functions, they are nested function calls)

# def fill(d, k):
#     if not k % 2:
#         if k/2 not in d:
#             fill(d,k/2)
#         print k
#         d[k] = 1 + d[k/2]
#     else:
#         o = 3*k+1
#         if o not in d:
#             fill(d,o)
#         print k
#         d[k] = 1 + d[o]
# d = {1:1}
# print d
# for i in range(2,1000000):
#     fill(d,i)
# highkey = 0
# highv = 0
# for (k,v) in d.iteritems():
#     if v > highv and k < 1000000:
#         highv = v
#         highkey = k
# print highkey
