from Simulation import *
from scipy.spatial import Delaunay

# 2d lubensky model 
# using framework

# setup for internal model is exactly the same 
# as in the 1d case
# only need to change cell locations and setup

# constants for model
Aa,As,Ah,Au = 0.65,0.5,1.5,2.2
na,ns,nh,nu = 4,4,8,8
ms,mh,mu = 4,4,6
Ts,Th,Tu = 4.0,101.0,2.0
Dh,Du = 200.0,0.16
S,H,U = 0.57,0.0088,4e-6
G,F = 0.8,0.6

# all cells have the same internal model
# 
# 4 species: 
IM = InternalModel()
IM.add_node('a','linear',params=[1.0])
IM.add_node('s','linear',params=[1.0/Ts])
IM.add_node('h','linear',params=[1.0/Th])
IM.add_node('u','linear',params=[1.0/Tu])

# internal interactions
# a -> a
IM.add_edge('a','a','hill_activ',params=[1.0,Aa,na])

# a -> s, s -> a
IM.add_edge('a','s','hill_activ',params=[1.0/Ts,As,ns])
IM.add_edge('s','a','hill_activ',params=[F,S,ms])

# a -> u
IM.add_edge('a','u','hill_activ',params=[1.0/Tu,Au,nu])

# a -> h
IM.add_edge('a','h','hill_activ',params=[1.0/Th,Ah,nh])
# h -> a
ha_edge = IM.add_edge('h','a','hill_activ',params=[G,H,mh])

# u -| (h -> a)
IM.add_edge('u',ha_edge,'hill_inactiv',is_mod=True,mod_type='mult',params=[1.0,1.0,U,mu])


# need to make some cells 
# we need to work a little harder in 2d
nrows = 10
ncolumns = 10
NCells = nrows * ncolumns

# first establish 'center' of each cell
# our indices i and j go like:
# i: left -> right
# j: bottom -> top
centers = np.empty((NCells,2)) # a center for each cell
for i in xrange(ncolumns):
    for j in xrange(nrows):
        c = np.array([float(i) + 0.5*(j % 2),(np.sqrt(3)/2)*float(j)])
        centers[i + j*ncolumns,:] = c

# can add noise to centers
# centers += np.random.normal(0.0,0.15,(NCells,2))        

# place each cell at these vertices
cells = [Cell(c) for c in centers]

# add these cells to the simulation
sim = Simulation()
for i in xrange(NCells):
    sim.add_cell(cells[i])

im_id = sim.add_internal_model(IM)

# set all cells to have the same internal model
sim.set_internal_model(range(NCells),im_id)

# need to figure out which cells are connect to which
connections = np.zeros((NCells,NCells)) > 0 # default boolean false array
tri = Delaunay(centers)
for k in xrange(NCells):
    # if index is found in triangle, add indices of points
    # in that triangle to connections
    matches = [t for t in tri.vertices if k in t]
    matches = reduce(lambda x,y: np.concatenate((x,y)),matches)
    # duplicates don't really matter since we're just setting them all to True
    connections[k,matches] = True

# add diffusion to h and u
sim.add_interaction('h','h','diffusion',connections,params=[Dh/Th])
sim.add_interaction('u','u','diffusion',connections,params=[Du/Tu])

# start with only first cell up
low_dict = {'a':0.0,'s':0.0,'h':0.0,'u':0.0}
high_dict = {'a':1.0+F,'s':1.0,'h':0.0,'u':0.0}
sim.set_initial_conditions(range(0,NCells),low_dict)
sim.set_initial_conditions([1,8],high_dict)

print 'starting simulation'
t = np.linspace(0,300,200)
cdata = sim.simulate(t)
print 'simulation done'


import matplotlib.pyplot as plt

# plot a at various times
tindex = 117
astatus = np.zeros(NCells)
for i in xrange(NCells):
    astatus[i] = cdata[i][tindex,0]
plt.scatter(centers[:,0],centers[:,1],
            c=astatus,s=50,marker='s',edgecolors='none')
plt.show()
