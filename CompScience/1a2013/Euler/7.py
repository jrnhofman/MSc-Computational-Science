#Crappy very slow method

# p = 2
# n = 5
# while p < 10001:
#     for x in range(3,n/3,2):
#         if n%x == 0:
#             n = n + 1
#             break
#         if x == n-2:
#             p = p + 1
#             print 'Prime',p,n
#             n = n + 2

# print 'The 10001st prime is',n

primelist = [2,3]
n = 5
p = 2
while p < 10001:
    for x in primelist:
        if n%x == 0:
            n = n + 2
            break
        elif x == primelist[-1]:
            primelist.append(n)
            p = p + 1
            n = n + 2
            break

print 'The 10001st prime is',n
