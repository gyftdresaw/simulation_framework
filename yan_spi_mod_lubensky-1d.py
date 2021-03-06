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

# new stuff
Ty,Tp,Tsp = 1.0,1.0,100

# all cells have the same internal model
# 
# 6 species: 
IM = InternalModel()
IM.add_node('a','linear',params=[1.0])
IM.add_node('s','linear',params=[1.0/Ts])
IM.add_node('h','linear',params=[1.0/Th])
IM.add_node('u','linear',params=[1.0/Tu])
# additional yan business
IM.add_node('y','linear',params=[1.0/Ty])
IM.add_node('p','linear',params=[1.0/Tp])
IM.add_node('sp','linear',params=[1.0/Tsp])

# internal interactions
# lubensky
# a -> a
aa_edge = IM.add_edge('a','a','hill_activ',params=[1.0,Aa,na])

# a -> s, s -> a
IM.add_edge('a','s','hill_activ',params=[1.0/Ts,As,ns])
IM.add_edge('s','a','hill_activ',params=[F,S,ms])

# a -> u
IM.add_edge('a','u','hill_activ',params=[1.0/Tu,Au,nu])

# a -> h
# no a -> h edge
# IM.add_edge('a','h','hill_activ',params=[1.0/Th,Ah,nh])
# h -> a
ha_edge = IM.add_edge('h','a','hill_activ',params=[G,H,mh])

# u -| (h -> a)
IM.add_edge('u',ha_edge,'hill_inactiv',is_mod=True,mod_type='mult',params=[1.0,1.0,U,mu])

# new edges
# replace a -> p
# ap_edge = IM.add_edge('a','p','hill_activ',params=[1.0/Tp,0.5,4])
# with a -> sp
IM.add_edge('a','sp','hill_activ',params=[4.0/Tsp,0.6,8])

# and sp -> p
IM.add_edge('sp','p','hill_activ',params=[6.0/Tp,1.0,1])

# p -> h
IM.add_edge('p','h','hill_activ',params=[0.5/Th,1.0,8])

# u -> y
uy_edge = IM.add_edge('u','y','hill_activ',params=[1.0/Ty,4e-6,6])

# completely speculative
IM.add_edge('h','sp','hill_activ',params=[2.0/Tsp,0.08,2])

# yan business interactions
# yan-pnt bistable switch:
#  yan -| pnt --> not yet
#  pnt -| yan
# IM.add_edge('y',ap_edge,'hill_inactiv',is_mod=True,mod_type='mult',params=[1.0,1.0,1.0,1])
IM.add_edge('p',uy_edge,'hill_inactiv',is_mod=True,mod_type='mult',params=[1.0,1.0,2.0,5])




# need to make some cells 
# the 1d case is easy:
# in our 'lattice', all the cells are distance 1 apart
NCells = 40
cells = [Cell([x]) for x in np.linspace(1,NCells,NCells)]

# add these cells to the simulation
sim = Simulation()
for i in xrange(NCells):
    sim.add_cell(cells[i])

im_id = sim.add_internal_model(IM)

# set all cells to have the same internal model
sim.set_internal_model(range(NCells),im_id)

# cells adjacent to one another are connected
# for diffusion we include main diagonal
# equivalent to 3 wide diagonal
diff_connections = (np.eye(NCells,k=-1) + np.eye(NCells,k=0) + np.eye(NCells,k=1)) > 0

# add diffusion to h and u
sim.add_interaction('h','h','diffusion',diff_connections,params=[Dh/Th])
sim.add_interaction('u','u','diffusion',diff_connections,params=[Du/Tu])

# spi diffusion
sim.add_interaction('sp','sp','diffusion',diff_connections,params=[0.5/Tsp])

'''
# no just lateral connections
# cells adjacent to one another are connected
# for activation, don't inlude main diagonal
hill_connections = (np.eye(NCells,k=-1) + np.eye(NCells,k=1)) > 0

# a -> pnt
sim.add_interaction('a','pnt','hill_activ',hill_connections,params=[1.0,1.5,4])
# a -> notch
sim.add_interaction('a','notch','hill_activ',hill_connections,params=[1.5,0.5,2])
'''
# start with only first cell up
low_dict = {'a':0.0,'s':0.0,'h':0.0,'u':0.0,'y':0.0,'p':0.0,'sp':0.0}
high_dict = {'a':1.0+F,'s':1.0,'h':0.0,'u':0.0,'y':0.0,'p':0.0,'sp':0.0}
sim.set_initial_conditions(range(1,NCells),low_dict)
sim.set_initial_conditions([0],high_dict)

print 'starting simulation'
t = np.linspace(0,200,500)
cdata = sim.simulate(t)
print 'simulation done'


import matplotlib.pyplot as plt
import matplotlib.cm as cm

x_coord = np.linspace(1,NCells,NCells)

# plot species at various times
def plot_species(species_name,times,x_coord,t):
    tindices = [np.abs(t-v).argmin() for v in times]
    astatus = np.zeros((NCells,len(tindices)))
    for i in xrange(NCells):
        for j in xrange(len(tindices)):
            astatus[i,j] = cdata[i][tindices[j],IM.get_node_id(species_name)]
    '''
    plt.figure()
    plt.scatter(x_coord,np.zeros(NCells),
                c=astatus,s=50,marker='s',edgecolors='none')
    plt.show()
    '''
    # plot along x axis
    plt.figure()
    colors = cm.Dark2(np.linspace(0, 1, len(tindices)))
    for j in xrange(len(tindices)):
        plt.scatter(x_coord,astatus[:,j],color=colors[j])
    plt.legend(['%.1f'% t for t in times],loc='best')
    plt.title(species_name)
    plt.show()

plot_species('a',np.linspace(0.0,150,8),x_coord,t)

def plot_cell(cid,(times),t):
    tstart,tend = [np.abs(t-v).argmin() for v in times]
    plt.figure()
    plt.plot(t[tstart:tend],cdata[cid][tstart:tend,:])
    plt.legend(IM.get_node_names())
    plt.show()
