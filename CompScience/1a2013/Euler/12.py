from math import sqrt,ceil

count = 2
tr_number = 1
add = 2

while True:
    for i in range(2,int(ceil(sqrt(tr_number)))):
        if tr_number%i == 0:
            count = count + 2
    print 'Triangle number',tr_number,', divisors',count
    if count >= 500:
        break
    tr_number = tr_number + add
    add = add + 1
    count = 2
