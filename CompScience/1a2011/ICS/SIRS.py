import random
from decimal import Decimal
from math import pi,exp,log
import cmath
from igraph import *

#Function making list of infected and recovered individuals and order parameter
def count_order():
    sigma = 0 #order parameter
    for i in range(N):
        
        #append if node is in infected state and add time step
        if g.vs[i]["time"] in range(t_inf,t_rec):
            infected.append(i)
	    if t>=t_cut:
            	phi = 2.*pi*(float(g.vs[i]["time"]-1.))/float(t_reset-1)
            	sigma += cmath.exp(1j*phi)
            g.vs[i]["time"] += 1

	#append if node is in recovered state and add time step
        elif g.vs[i]["time"] in range(t_rec,t_reset):
	    recover.append(i)
	    if t>=t_cut:
		phi = 2.*pi*(float(g.vs[i]["time"]-1.))/float(t_reset-1)
        	sigma += cmath.exp(1j*phi)
	    g.vs[i]["time"] += 1
        
    #return S,I and R nodes and sigma
    if (len(infected)+len(recover))!=0:
        return N-len(infected)-len(recover),len(infected),len(recover),abs(sigma/(len(infected)+len(recover)))
    else:
        return N-len(infected)-len(recover),len(infected),len(recover),0

#Function infect neighbours of infected nodes
def infect():
    for i in range(len(infected)):
        
        #calculate neighbors of node
    	nb = g.neighbors(infected[i])
   	for j in range(len(nb)):
	    if g.vs[nb[j]]["time"] == 0:
            	r = random.random()

                #infect neighbor with probability p_inf
                if r < p_inf: 
                    g.vs[nb[j]]["time"] = t_inf
    return None

#Parameters
p_inf = 0.1 #probability of infection
t_max = 3000 #maximum run time
t_cut = 1000 #cut off time for calculating sigma
N = 1000 #size of graph
K = 5 #number of connections per node for p = 0
real = 50 #number of graph realizations per p-value

#Time parameters for nodes
t_susc = 0 
t_inf = 1 #infected steps 9-1 = 8
t_rec = 9 #recovered steps 17-9 = 8
t_reset =17 #total cycle time

#Open file for writing
f1 = file('/home/jhofman/Desktop/Project/sigma_test_test.txt','a')

#Construction of list for p-values
plist1 = [i*0.005 for i in range(1,20)]
plist2 = [2**i for i in range(-15,1)]
plist3 = [i/20. for i in range(0,21)]
plist2.extend(plist3)
plist1.extend(plist2)

#Loop over all p
for p in plist1:

    sigma_mean = 0
    
    #loop over number of realizations
    for m in range(real):
        
        sigma = 0 

        #creates graph with Watts-Strogatz algorithm and remove loops
        g = Graph.Watts_Strogatz(1,N,K,p)
        g = Graph.simplify(g) 
        
        #set all nodes to susceptible state except for the random seed
        for i in range(N):
            g.vs[i]["time"] = 0
        seed = random.randint(0,N-1)
        g.vs[seed]["time"] = 1

        #loop over t
        t = 0
        
        while t < t_max:
        
            infected = []
	    recover = []
	    (susc,inf,rec,sig) = count_order() 
	    sigma += sig #sum up all the order parameters for all time steps
            infect() 
            
            #put some recovered nodes back in susceptible state
            for i in range(len(recover)):
		g.vs[recover[i]]["time"] = (g.vs[recover[i]]["time"]%t_reset) 
            
            t += 1

        #calculates the order parameter over time t_max-t_cut
        sigma_mean += sigma/float(t_max-t_cut) 
        print "p = ",p," m = ",m," sigma = ",sigma_mean

    #save the order parameter as a mean over all realizations
    f1.write(str(Decimal(p))+' '+str(sigma_mean/float(real))+'\n') 
    print "p = ",p," sigma_mean_mean = ",sigma_mean/float(real)
    
f1.close()
