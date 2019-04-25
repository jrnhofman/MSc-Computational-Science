from random import expovariate,randint,random
from collections import deque
from numpy import mean,std
from math import sqrt,floor
#import matplotlib.pyplot as plt
#import matplotlib.mlab as lab

"""
Variables
"""
maxTime = 1000000

"""
Fixed parameters
"""
maxCapacity = 6
boardTime = 2 #also getting off the elevator
doorOpening = 2
doorClosing = 2
accTime = 1
betwFloorTime = 3
nFloors = 4
c = 0.9

"""
Initialising queues and arrays
"""
meanTravelTime = list()
meanTravelTimeError = list()
meanSampleDiff = list()
meanSampleDiffError = list()

"""
Calculates next arrival time
"""
def arrival(time,MTBA):
    random = expovariate(1./MTBA)
    time += random
    return time

"""
Calculates departure time of elevator
"""
def departure(time):
    if(time==nextArrival): #if someone just arrived at elevator
        if(nextDeparture < 10**7): #if departure already scheduled, add board time
            time = nextDeparture + boardTime
        else: #else open door + board + close door
            time += doorOpening + boardTime + doorClosing
    elif(time==nextReturn and len(queue)>0): #if the elevator just returned and people are in queue
        time += doorOpening + min(boardTime*maxCapacity,boardTime*len(queue)) + doorClosing
    else: #do not schedule departure
        time = 10**7
    return time

"""
Calculates return time of elevator
"""
def returning(time,capacity):
    departTime = time
    capacity = deque(sorted(list(capacity),key = lambda dest: dest[2])) #sort according to floor destination
    time += capacity[0][2] * betwFloorTime + accTime * 2 #getting to floor of first person
    time += doorOpening + boardTime #open door and exit first person
    if(measure) : travelTime.append(time - departTime)
    prevLoc = capacity[0][2]
    capacity.popleft()
    while(len(capacity)>0): #do for the rest of the people in elevator
        if capacity[0][2] == prevLoc: #exit when on same floor
            time += boardTime
        else: #close door, switch floor, open door and exit
            time += doorClosing
            time += (capacity[0][2] - prevLoc) * betwFloorTime + accTime * 2
            time += doorOpening
            time += boardTime
        if(measure) : travelTime.append(time - departTime)
        prevLoc = capacity[0][2]
        capacity.popleft()
    time += doorClosing 
    time += prevLoc * betwFloorTime + accTime * 2 #head back to ground floor
    return time,capacity


"""
Loop over MTBA list
"""
MTBAList = [x+5 for x in range(30)]
for m in MTBAList:

    queueTime = list()
    travelTime = list()
    journeyTime = list()
    sampleDiff = list()
    queue = deque()
    capacity = deque()

    MTBA = m

    time = 0.
    nextArrival = arrival(time,MTBA) #set first arrival
    nextDeparture = 10**7 #this means none scheduled
    nextReturn = 10**7 #this means none scheduled
    nDepartures, nArrivals = 0,0

    measure = False

    """
    Run simulation
    """
    while(time<maxTime):
        """
        Measure every other 1000 timesteps
        and save average journey time and 
        sample difference after every
        measurement period
        """
        if(floor(time/1000)%2==0 and measure == True): 
            measure = False
            totalTime = mean(travelTime) + mean(queueTime)
            journeyTime.append(totalTime)
            if(time>3000): sampleDiff.append(oldMean-totalTime)
            oldMean = totalTime
        if(floor(time/1000)%2==1 and measure == False):
            measure = True
            queueTime = list()
            travelTime = list()
            

        """
        Customer arrives
        """
        if(nextArrival <= nextDeparture and nextArrival <= nextReturn):
            time = nextArrival
            dest = randint(1,4)
            #Case 2, people may skip with dest floor 1 
            #if(random()>(1-c**len(queue)) or dest!=1):#passing probability
            
            #Case 3, people may always skip
            if(random()>(1-c**len(queue))/dest):
            
            #Case 1, don't skip    
            #if(0==0):    
                queue.append([nArrivals,time,dest])
                nArrivals += 1
                #Only alter departure time if there is room and elevator is present
                if(len(queue)<=maxCapacity and nextReturn == 10**7):
                    nextDeparture = departure(time) 
            #print 'Arrival, {0:10f} {1:4d} {2:4d} {3:1d}'.format(time,nArrivals,len(queue),queue[-1][2])
            nextArrival = arrival(time,MTBA)

        """
        Elevator takes off
        """
        if(nextDeparture < nextArrival and nextDeparture <= nextReturn):
            time = nextDeparture
            for i in range(min(maxCapacity,len(queue))):
                if(measure) : queueTime.append(time - queue[0][1])
                capacity.append(queue[0])
                queue.popleft()
                nDepartures += 1
            #print 'Departu, {0:10f} {1:4d} {2:4d}'.format(time,nDepartures,len(queue))
            (nextReturn,capacity) = returning(time,capacity)
            nextDeparture = 10**7
        
        """
        Elevator returns
        """
        if(nextReturn < nextArrival and nextReturn < nextDeparture):
            time = nextReturn
            #print 'Returns, {0:10f}'.format(time)
            nextDeparture = departure(time)
            nextReturn = 10**7
    
    """
    Arrays for the mean travel time and difference between
    subsequent measurement period, one entry per MTBA value
    """
    meanTravelTime.append(mean(journeyTime))
    meanTravelTimeError.append(1.96*std(journeyTime)/sqrt(maxTime/2000))
    meanSampleDiff.append(mean(sampleDiff))
    meanSampleDiffError.append(1.96*std(sampleDiff)/sqrt(maxTime/2000-1))
    print "\n"
    print "Mean time between arrivals = ",m
    print "Overall mean journey time = ",mean(journeyTime)
    print "Error = ",1.96*std(journeyTime)/sqrt(maxTime/2000)
    print "Average difference between samples = ",mean(sampleDiff)
    print "Error = ",1.96*std(sampleDiff)/sqrt(maxTime/2000-1)


"""
Plotting and saving
"""
#lab.save("ElevatorSamplePassMult.dat",(MTBAList,meanSampleDiff,meanSampleDiffError))
#lab.save("ElevatorPassMult.dat",(MTBAList,meanTravelTime,meanTravelTimeError))
#plt.figure(1)
#plt.errorbar(MTBAList,meanSampleDiff,meanSampleDiffError,fmt='b.')
#plt.xlabel("Mean time between arrivals")
#plt.ylabel("Average difference between samples")
#plt.figure(2)
#del MTBAList[0:3]
#del meanTravelTime[0:3]
#del meanTravelTimeError[0:3]
#plt.errorbar(MTBAList,meanTravelTime,meanTravelTimeError,fmt='rs')
#plt.xlabel("Mean time between arrivals")
#plt.ylabel("Average journey time")
#plt.show()







