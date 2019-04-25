shift = [3,0,3,2,3,2,3,3,2,3,2,3]
shiftleap = [3,1,3,2,3,2,3,3,2,3,2,3]

year = 1
startday = 1+sum(shift) #1900 is not a leap year
sum = 0

for year in range(1,101):
    for i in range(len(shift)):
        if year%4!=0:
            startday += shift[i]
        else:
            startday += shiftleap[i]
        if startday%7==0:
            sum += 1
print 'The number of Sundays on the first of the month from 1 Jan 1901 to 31 Dec 2000 is',sum
