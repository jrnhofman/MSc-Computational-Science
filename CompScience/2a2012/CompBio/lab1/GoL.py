import pylab as PL
import random as RD
import scipy as SP

RD.seed()

# set parameters
width = 50
height = 50
initProb = 0.2
maxTime = 1000

# initialize configuration
config = SP.zeros([height, width])
config[width/2,height/2] = 1
config[width/2-1,height/2+1] = 1
config[width/2,height/2-1] = 1
config[width/2-1,height/2] = 1
config[width/2+1,height/2] = 1

# visualize configuration
PL.ion()
fig = PL.figure()
PL.hold(False)
PL.axis('image')

def visualize(data, time):
    PL.cla()
    PL.pcolor(data, vmin = 0, vmax = 1, cmap = PL.cm.binary)
    PL.title('t = ' + str(time))
    fig.canvas.manager.window.update()

visualize(config, 0)

# main loop
nextConfig = SP.zeros([height, width])
for t in xrange(maxTime):
    # simulate Game of Life rule
    for x in xrange(width):
        for y in xrange(height):
            state = config[y, x]
            numberOfAlive = 0
            for dx in xrange(-1, 2):
                for dy in xrange(-1, 2):
                    numberOfAlive += config[(y+dy)%height, (x+dx)%width]
            if state == 0 and numberOfAlive == 3:
                state = 1
            elif state == 1 and (numberOfAlive < 2 or numberOfAlive > 3):
                state = 0
            nextConfig[y, x] = state
    # swap names of configs
    config, nextConfig = nextConfig, config
    visualize(config, t + 1)
fig.canvas.manager.window.wait_window()
