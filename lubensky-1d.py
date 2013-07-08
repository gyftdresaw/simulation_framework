from Simulation import *

# 1d lubensky model 
# using framework

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
# the 1d case is easy:
# in our 'lattice', all the cells are distance 1 apart
NCells = 50
cells = [Cell([x]) for x in np.linspace(1,NCells,NCells)]

# add these cells to the simulation
sim = Simulation()
for i in xrange(NCells):
    sim.add_cell(cells[i])

im_id = sim.add_internal_model(IM)

# set all cells to have the same internal model
sim.set_internal_model(range(NCells),im_id)

# cells adjacent to one another are connected
# equivalent to 3 wide diagonal
connections = (np.eye(NCells,k=-1) + np.eye(NCells,k=0) + np.eye(NCells,k=1)) > 0

# add diffusion to h and u
sim.add_interaction('h','h','diffusion',connections,params=[Dh/Th])
sim.add_interaction('u','u','diffusion',connections,params=[Du/Tu])

# start with only first cell up
low_dict = {'a':0.0,'s':0.0,'h':0.0,'u':0.0}
high_dict = {'a':1.0+F,'s':1.0,'h':0.0,'u':0.0}
sim.set_initial_conditions(range(1,NCells),low_dict)
sim.set_initial_conditions([0],high_dict)

print 'starting simulation'
t = np.linspace(0,300,1000)
cdata = sim.simulate(t)
print 'simulation done'


import matplotlib.pyplot as plt

# plot a at various times
tindex = 500
astatus = np.zeros(NCells)
for i in xrange(NCells):
    astatus[i] = cdata[i][tindex,0]
plt.scatter(np.linspace(1,NCells,NCells),np.zeros(NCells),
            c=astatus,s=50,marker='s',edgecolors='none')
