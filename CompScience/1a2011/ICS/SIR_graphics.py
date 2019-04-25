from igraph import *
import random

#Video options (standard disabled)
from video import MEncoderVideoEncoder
encoder = MEncoderVideoEncoder()
encoder.add("Sample epidemic spreading")
encoder.fps = 2

#Function for setting state and counting states
def count():
    susc,inf,rec = 0,0,0
    for k in range(0,length):
        if g.vs[k]["time"] in range(t_susc,t_inf):
            g.vs[k]["state"] = 'susc'
            susc += 1
        elif g.vs[k]["time"] in range(t_inf,t_rec):
            inf += 1
            g.vs[k]["state"] = 'inf'
        else:
            rec += 1
            g.vs[k]["state"] = 'rec'
    return susc, inf, rec

#Parameters
p = 0.5
p_inf = 0.1
t_max = 1000

t_susc = 0
t_inf = 1
t_rec = 5


g = Graph.Watts_Strogatz(1,100,3,p)
g = Graph.simplify(g)
layout = g.layout("fruchterman_reingold")
colour_dict = {"susc" : "white","inf" : "red","rec" : "green"}
length = len(g.vs)
for i in range(0,length):
    g.vs[i]["time"] = 0

seed = random.randint(0,length-1)
g.vs[seed]["time"] = 1

#Iteration step
t = 0
while t < t_max:
    for i in range(0,length):
        if g.vs[i]["time"] in range(t_inf,t_rec):
            g.vs[i]["time"] += 1
            nb = g.neighbors(i)
            for j in range(0,len(nb)):
                r = random.random()
                if r < p_inf and g.vs[nb[j]]["time"] == 0:
                    g.vs[nb[j]]["time"] = t_inf
    (susc,inf,rec) = count()
    #print "t = ",t," suscept = ",susc," inf = ",inf," recov = ",rec
    g.vs["color"] = [colour_dict[state] for state in g.vs["state"]]
    q = Plot(bbox=(500,400))
    q.add(g,layout = layout,vertex_label=None,legend=True)
    encoder.add(q)
    t += 1

#Plot final state
s = Plot(bbox=(500,400))
layout = g.layout("fruchterman_reingold")
s.add(g,layout = layout,vertex_label=None)
#s.add('hoi')
s.show()

#Save video (standard disabled)
encoder.save("sample_SIR.avi")
