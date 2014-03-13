from Simulation import *

# 1d lubensky model 
# using framework
# DEMONSTRATING SHORTENED LATTICE SPACING

# constants for model
Aa,As,Ah,Au = 0.60,0.5,1.5,2.2 # reducing A self activation threshold to pull up faster
na,ns,nh,nu = 4,4,8,8
ms,mh,mu = 4,4,6
Ts,Th,Tu = 4.0,101.0,2.0
Dh,Du = 200.0,0.04 # Reducing diffusion of inhibitor u to decrease lattice
# reducing Du to 0.02 is TOO short, start to see Ato instability (twinning)
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
NCells = 30
cells = [Cell([x]) for x in np.linspace(1,NCells,NCells)]

# add these cells to the simulation
sim = Simulation()
for i in xrange(NCells):
    sim.add_cell(cells[i])

im_id = sim.add_internal_model(IM)

# set all cells to have the same internal model
sim.set_internal_model(range(NCells),im_id)

# set boundary condition
sim.set_boundary_conditions([0],'ref_on')

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
t = np.linspace(0,150,300)
cdata = sim.simulate(t)
print 'simulation done'


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

x_coord = np.linspace(1,NCells,NCells)

# plot species at various times
def plot_species(species_name,times):
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
    plt.legend(['%.1f'% time for time in times],loc='best')
    plt.title(species_name,fontsize=24)
    # plt.show()

def plot_all_species(times):
    for n in IM.node_names.keys():
        plot_species(n,times)

plot_species('a',np.linspace(0.0,150,8))

def plot_cell(cid,times):
    tstart,tend = [np.abs(t-v).argmin() for v in times]
    plt.figure()
    plt.plot(t[tstart:tend],cdata[cid][tstart:tend,:])
    plt.legend(IM.get_node_names())
    # plt.show()

# differentiate cell types
def get_r8s(time):
    tindex = np.abs(t-time).argmin()
    return [i for i in xrange(NCells) if cdata[i][tindex,IM.get_node_id('u')] > 0.010]

def get_r25s(time):
    r25_indices = [ [r8-1,r8+1] for r8 in get_r8s(time)]
    # not the most readable just flattens the list
    return [c for pair in r25_indices for c in pair if c >= 0 and c < NCells]

def get_progenitors(time):
    r8s = get_r8s(time)
    r25s = get_r25s(time)
    return [c for c in xrange(NCells) if not c in r8s and not c in r25s]


## plot based on these categorizations

def plot_cell_species(species_name,time):
    tindex = np.abs(t-time).argmin()
    astatus = np.zeros((NCells,1))
    for i in xrange(NCells):
        astatus[i] = cdata[i][tindex,IM.get_node_id(species_name)]
    r8s = get_r8s(time)
    r25s = get_r25s(time)
    progenitors = get_progenitors(time)
    cell_types = [r8s,r25s,progenitors]
    # plot along x axis
    plt.figure()
    gs = gridspec.GridSpec(2,1,height_ratios=[19,1])
    # first plot values 
    ax1 = plt.subplot(gs[0])
    ax1.xaxis.grid(True,linewidth=0.1,linestyle='-',color='0.4');
    colors = cm.Set1(np.linspace(0, 1, 10))
    for j in xrange(len(cell_types)):
        plt.scatter(x_coord[cell_types[j]],astatus[cell_types[j]],color=colors[j],s=50)
    plt.legend(['R8','R2/5','MP'],loc='best')
    for j in xrange(len(cell_types)):
        plt.plot(x_coord[cell_types[j]],astatus[cell_types[j]],'--',color=colors[j])
    plt.title(species_name,fontsize=24)
    plt.xlabel('Cell',labelpad=0)
    # now plot stripes
    ax2 = plt.subplot(gs[1],sharex=ax1)
    xcorners = np.tile(np.linspace(0,NCells,NCells+1)+0.5,(2,1))
    ycorners = np.tile(np.array([[0],[1]]),(1,NCells+1))
    plt.pcolor(xcorners,ycorners,astatus.T,edgecolors='k')
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    # plt.show()

def plot_all_cell_species(time):
    for n in IM.node_names.keys():
        plot_cell_species(n,time)
