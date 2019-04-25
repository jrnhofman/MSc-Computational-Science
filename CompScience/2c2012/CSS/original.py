from numpy import *
from random import *
import matplotlib.pyplot as plt
from matplotlib import colors
import os, sys
import time

# 0 empty side
# 1 inactive agent
# 2 active agent
# 3 jailed agent
# 4 cop

##############
# Parameters #
##############
MAKE_AWESOME_MOVIE = False

gridSize = 40
agentDensity = 0.7
copDensity = 0.04
legitimacy = 0.8
agentVision = 3
copVision = 3
maxJail = 15
threshold = 0.1
k = 2.3
maxIter = 1000
visionLength = 2*agentVision+1

class Agent:

    def __init__(self,legitimacy,vision,grid,ids):
        self.hardship = uniform(0,1)
        self.legitimacy = legitimacy
        self.grievance = self.hardship * (1-self.legitimacy)
        self.risk_aversion = uniform(0,1)
        self.vision = vision
        self.active = 1
        self.id = ids
        self.maxJail = -1
        self.jail = 0

        #Place agent randomly
        position = (randint(0,gridSize-1),randint(0,gridSize-1))
        while (grid[position] != 0) and (grid[position] != 3):
            position = (randint(0,gridSize-1),randint(0,gridSize-1))
        self.position = position
        grid[self.position] = self.active

    def act(self,vision):

        #If not jailed:
        if self.maxJail == -1:
            xx = self.position[0]
            yy = self.position[1]

            #Move to random unoccupied position
            grid[self.position] = 0
            position = ((xx+randint(-vision,vision))%gridSize,(yy+randint(-vision,vision))%gridSize)
            while (grid[position] != 0) and (grid[position] != 3):
                position = ((xx+randint(-vision,vision))%gridSize,(yy+randint(-vision,vision))%gridSize)
            self.position = position
            grid[self.position] = self.active

            #Compute active to cop ratio in vision
            arr=roll(roll(grid,shift=-self.position[0]+agentVision,axis=0),shift=-self.position[1]+agentVision,axis=1)[:visionLength,:visionLength].flatten()

            if arr.max()==4:

                cops = bincount(arr)[4]
                active = bincount(arr)[2]

                if self.active == 1:
                    active += 1

                #Arrest probability and net risk
                P = 1 - exp(-k*floor(cops/active))
                self.net_risk = self.risk_aversion * P
            else:
                self.net_risk = self.risk_aversion * 0.9

            #Become active or inactive
            if (self.grievance - self.net_risk) > threshold:
                self.active = 2
                grid[self.position] = 2
            else:
                self.active = 1
                grid[self.position] = 1

        #If jailed:
        elif self.jail < self.maxJail-1:
            self.jail += 1

        #If at last jail-turn:
        else:
            self.maxJail = -1
            self.active = 1
            grid[self.position] = 1

class Cop:

    def __init__(self,vision,grid):
        self.vision = vision

        position = (randint(0,gridSize-1),randint(0,gridSize-1))
        while (grid[position] != 0) and (grid[position] != 3):
            position = (randint(0,gridSize-1),randint(0,gridSize-1))
        self.position = position
        grid[self.position] = 4


    def act(self,vision):

        #Move cop to random non-occupied grid point in vision
        xx = self.position[0]
        yy = self.position[1]

        grid[self.position] = 0
        position = ((xx+randint(-vision,vision))%gridSize,(yy+randint(-vision,vision))%gridSize)
        while (grid[position] != 0) and (grid[position] != 3):
            position = ((xx+randint(-vision,vision))%gridSize,(yy+randint(-vision,vision))%gridSize)
        self.position = position
        grid[self.position] = 4


        arr=roll(roll(grid,shift=-self.position[0]+agentVision,axis=0),shift=-self.position[1]+agentVision,axis=1)[:visionLength,:visionLength]
        if len(where(arr==2)[0])>0:

            #Make list of active agents in vision
            xx = self.position[0]
            yy = self.position[1]
            active = []
            for x in range(-vision,vision+1):
                hor = (xx+x)%gridSize
                for y in range(-vision,vision+1):
                    ver = (yy+y)%gridSize
                    if grid[hor,ver] == 2:
                        active.append((hor,ver))

            #If there are any active agents, arrest one and move there:
            grid[self.position] = 0
            pos = active[randint(0,len(active)-1)]
            b = [x for x in Agents if x.position==pos]
            b[0].active = 3
            b[0].jail = 0
            b[0].maxJail = randint(1,maxJail)
            self.position = pos
            grid[pos] = 4


#Keep track of jailed, inactive and active per iteration
def count(iters):
    active,inactive,jailed = 0,0,0
    for x in Agents:
        if x.active == 1:
            inactive += 1
        elif x.active == 2:
            active += 1
    ac.append([iters,active])
    ina.append([iters,inactive])
    ja.append([iters,agentDensity*gridSize*gridSize-active-inactive])


#Make a plot of the current grid
def visualize():
    gridv = zeros((gridSize,gridSize),dtype=int)
    for x in Agents:
        gridv[x.position] = x.active
    for x in Cops:
        gridv[x.position] = 4
    return gridv


#Relevant statistics
def stat(fileac,fileinac,fileja):
    ac = loadtxt('Figures/'+fileac+'.dat')
    inac = loadtxt('Figures/'+fileinac+'.dat')
    ja = loadtxt('Figures/'+fileja+'.dat')

    fig2 = plt.figure(0)
    plt.plot(delete(ac,0,1)[:1000],'red',label='Active Agents')
    plt.plot(delete(inac,0,1)[:1000],'green',label='Inactive Agents')
    plt.plot(delete(ja,0,1)[:1000],'gray',label='Jailed Agents')
    plt.ylabel('Count')
    plt.xlabel('Time')
    plt.legend()
    fig2.show()

    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)
    peak = [x for x in ac if x[1]>50]
    count = 1
    diff,total,length = [],[],[]
    for x in range(len(peak)-1):
        if (peak[x+1][0] - peak[x][0]) > 1:
            diff.append(peak[x+1][0]-peak[x][0])
            length.append(count)
            count = 1
        else:
            count += 1

    ax1.hist(diff,bins=50,color='green')
    ax1.set_ylabel('Count')
    ax1.set_xlabel('Waiting Time')
    ax2 = fig.add_subplot(212)
    summ = 0
    for x in range(len(peak)-1):
        if (peak[x+1][0] - peak[x][0]) > 1:
            summ += peak[x][1]
            total.append(summ)
            summ = 0
        else:
            summ += peak[x][1]
    ax2.hist(total,bins=50,color='green')
    ax2.set_ylabel('Count')
    ax2.set_xlabel('Total Activition')
    fig.show()

    bins = [i for i in range(21)]
    fig3 = plt.figure(2)
    plt.hist(length,bins=bins,color='green')
    plt.ylabel('Count')
    plt.xlabel('Outburst Duration')
    fig3.show()

    print 'Waiting time: \t \t','Mean: ',mean(diff,axis=0),'\t Std: ',std(diff,axis=0),'\t Median: ',median(diff,axis=0)
    print 'Total Activation:\t','Mean: ',mean(total,axis=0),'\t Std: ',std(total,axis=0),'\t Median: ',median(total,axis=0)
    print 'Outburst Duration:\t','Mean: ',mean(length,axis=0),'\t Std: ',std(length,axis=0),'\t Median: ',median(length,axis=0)
    print len(total), len(diff)

    return diff, total, peak


#Initialize all the agents and cops
grid = zeros((gridSize,gridSize),dtype=int)
Agents = []
for i in range(int(agentDensity*gridSize*gridSize)):
    Agents.append(Agent(legitimacy,agentVision,grid,i))
Cops = []
for i in range(int(copDensity*gridSize*gridSize)):
    Cops.append(Cop(copVision,grid))


if MAKE_AWESOME_MOVIE:
    # Some figure initializations
    cmap = colors.ListedColormap(['white','green','red','blue','black'])
    norm = colors.BoundaryNorm([-0.5 + i for i in range(6)], cmap.N)
    plt.figure()

    plt.matshow(visualize(),cmap=cmap,norm=norm)
    fname = 'Figures/_tmp%03d.png'%0
    plt.savefig(fname)

########
# Main #
########

start = time.time()
ac,ina,ja = [],[],[]
count(0)

#Main loop
for iters in range(maxIter):
    if iters%100==0:
        print 'Iteration: ',iters

    #Create shuffled lists of all objects
    x = range(len(Agents+Cops))
    shuffle(x)

    #Go over list and agents one by one
    for i in x:
        if i < len(Agents):
            Agents[i].act(agentVision)
        else:
            Cops[i-len(Agents)].act(copVision)
    count(iters+1)

    if MAKE_AWESOME_MOVIE:

        # Saving figures and such
        plt.clf()
        plt.matshow(visualize(),cmap=cmap,norm=norm)

        fname = 'Figures/_tmp%03d.png'%(iters+1)
        print 'Saving frame',fname
        plt.savefig(fname)


if MAKE_AWESOME_MOVIE:
    # Making a movie
    print 'Making movie animation.mpg - this make take a while'
    os.system("mencoder 'mf://Figures/_tmp*.png' -mf type=png:fps=3 -ovc lavc -oac copy -o original_+"+str(agentVision)+".mpg")

print 'Time per iteration: ',(time.time() - start)/(maxIter+1)

#Save some data
# savetxt('tmp/ac'+str(agentVision)+'.dat',ac)
# savetxt('tmp/ina'+str(agentVision)+'.dat',ina)
# savetxt('tmp/ja'+str(agentVision)+'.dat',ja)

#plt.show(plt.plot(delete(loadtxt('Figures/ac7.dat'),0,1)))

