from igraph import *
import random
from decimal import Decimal

#Open files for writing
f1 = file('/home/jhofman/Desktop/Project/rhalve.txt','a')
f2 = file('/home/jhofman/Desktop/Project/rtotal.txt','a')

#Function for setting state and counting states
def count():
    for i in range(0,N):
        
        #Append to infected list and add time step
        if g.vs[i]["time"] in range(t_inf,t_rec):
            infected.append(i)
            g.vs[i]["time"] += 1

        #Append to recovered list and add time step
        elif g.vs[i]["time"] == t_rec:
            recover.append(i)

    #Return S,I and R counts and the infected+recovered portion
    return 1000-(len(infected)+len(recover)),len(infected), len(recover),(len(infected)+len(recover))/1000.

#Function doing the infection step 
def infect():
    for i in range(len(infected)):
        
        #Calculate neighbors of infected node
        nb = g.neighbors(infected[i])
        for j in range(len(nb)):
            if g.vs[nb[j]]["time"] == 0:
                r = random.random()
                
                #Infect with probability p_inf
                if r <= p_inf:
                    g.vs[nb[j]]["time"] = t_inf
    return None

#Parameters
p_inf = 1.0 #Probability of infection
t_max = 1000 #Maximum time of run
n = 100  #number of realizations
N = 1000 #number of nodes

t_susc = 0
t_inf = 1 #infection time 5-1 = 4
t_rec = 5

plist = [2**i for i in range(-15,1)] #list of p-values
for p in plist:
    r_total_sum,r_half_sum = 0.,0.
    
    #loop over n realizations
    for m in range(n):
        t,check,r_half = 0,0,0
        
        #Initialize graph and simplify (remove double loops)
        g = Graph.Watts_Strogatz(1,1000,5,p)
        g = Graph.simplify(g)
        
        #Set all nodes to susceptible state except random seed
        for i in range(N):
            g.vs[i]["time"] = 0
        seed = random.randint(0,N-1)
        g.vs[seed]["time"] = 1

        #Iteration step
        while r_half < 1.0:
            infected = []
            recover = []
           
            (susc,inf,rec,r_half) = count()
            infect()
            
            #record time to infect halve the network
            if r_half > 0.5 and check == 0:
                half_time = t
                r_half_sum += t
                check = 1
            
            #record time to infect the whole network    
            if r_half == 1.:
                r_total_sum += t
            
            #if time reaches t_max while r_half<1 some nodes have not been infected, so generate new graph
            if t==t_max:
                m -= 1
                r_half_sum -= half_time
                break

            t += 1
        
    #write averaged half and whole infection rates to file    
    r_half_sum /= n
    r_total_sum /= n
    print "t_rec = ",t_rec-1," p = ",p," r_half_sum = ",r_half_sum
    print "t_rec = ",t_rec-1," p = ",p," r_total_sum = ",r_total_sum
    f1.write(str(t_rec-1)+' '+str(Decimal(p))+' '+str(r_half_sum)+'\n')
    f2.write(str(t_rec-1)+' '+str(Decimal(p))+' '+str(r_total_sum)+'\n')

f1.close()
f2.close()
