from igraph import *
import random
from decimal import Decimal
import numpy
from numpy import std

#Open file for writing
f = open('/home/jhofman/Desktop/Project/test.txt','a')

#Function for calculating global clustering coefficient
def connectivity(g):

    #Gives a (symmetric) matrix of connections between nodes
    matrix = g.get_adjacency()
    conn = 0
    for i in range(len(g.vs)):
        p = 0
        
        #List the neighbors of node i
        list = g.neighbors(i)
        if len(list)>1:
            
            #Maximum possible number of connections
            max_p = len(list)*(len(list)-1.)/2.

            #Number of connections among neighbors using the adjacency matrix
            for j in range(len(list)):
                for k in range(j+1,len(list)):
                    if matrix[list[j]][list[k]]==1:
                        p +=1
            
            #add ratio to average connectivity
            conn += p/max_p

        else: conn += 1
    return conn/len(g.vs)

#Parameters
n = 100 #number of realizations
c0,l0 = 0,0
stdd = 0

#loop over realizations
for m in range(n):

    #Initialize graph with rewiring probability 0 and compute transitivity (=connectivity) and average path length
    g0 = Graph.Watts_Strogatz(1,1000,5,0)
    c0 += connectivity(g0) 
    l0 += g0.average_path_length()
    trans = GraphBase.transitivity_local_undirected(g0,g0.vs)

    #standarddeviation of transitivity
    stdd += numpy.std(trans)

#Print and write averages over realizations
print "p = ",0," Relative average path length = ",l0/l0," Relative connectivity =",c0/c0," std = ",stdd/float(n) 
f.write('p = '+str(0)+' l(p) = '+str(l0/l0)+' c(p) = '+str(c0/c0)+' '+str(Decimal(stdd/float(n)))+'\n')    

#Loop over values in p-lists
plist = [2**i for i in range(-15,1)]
p2list = range(0.07,0.25,0.01)
p2list.append(range(0.3,1,0.1))
p2list.extend(plist)
for p in p2list:
    
    stdd = 0
    c,l = 0,0	
    m = 0
    
    #loop ove realizations
    while m < n:
        
        #Initialize graph and calculate properties
        g = Graph.Watts_Strogatz(1,1000,5,p)
        trans = GraphBase.transitivity_local_undirected(g,g.vs)
	
        #Make sure graph is connected
        if not math.isnan(numpy.std(trans)):
            c += connectivity(g)
            l += g.average_path_length()
            stdd += numpy.std(trans)
	    m += 1	
    
    #Average properties over realizations and normalize over p = 0 properties
    print "p = ",Decimal(p)," Relative average path length = ",l/float(l0)," Relative connectivity =",c/float(c0)," std = ",stdd/float(n)
    f.write('p = '+str(Decimal(p))+' l(p) = '+str(l/float(l0))+' c(p) = '+str(c/float(c0))+' '+str(Decimal(stdd/float(n)))+'\n')

f.close()

    
