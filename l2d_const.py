from Simulation import *
from scipy.spatial import Delaunay

# 2d lubensky model 
# using framework
# 
# first row of cells is given a constant source of h to propagate R8's

# setup for internal model is exactly the same 
# as in the 1d case
# only need to change cell locations and setup

# constants for model
Aa,As,Ah,Au = 0.65,0.5,1.5,2.2 # decrease Ah threshold from 1.5
na,ns,nh,nu = 4,4,8,8
ms,mh,mu = 4,4,6
Ts,Th,Tu = 4.0,101.0,2.0
Dh,Du = 200.0,0.16 # reducing Du by a little from 0.16
S,H,U = 0.57,0.0088,4e-6 # decrease activation of H from 0.0088
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

## h source cells have an internal model without h connections/perturbations ##
IM_source = InternalModel()
IM_source.add_node('a','linear',params=[1.0])
IM_source.add_node('s','linear',params=[1.0/Ts])
IM_source.add_node('h') # MODIFIED no h degradation
IM_source.add_node('u','linear',params=[1.0/Tu])

# internal interactions
# a -> a
IM_source.add_edge('a','a','hill_activ',params=[1.0,Aa,na])

# a -> s, s -> a
IM_source.add_edge('a','s','hill_activ',params=[1.0/Ts,As,ns])
IM_source.add_edge('s','a','hill_activ',params=[F,S,ms])

# a -> u
IM_source.add_edge('a','u','hill_activ',params=[1.0/Tu,Au,nu])
# and no h related interactions

# need to make some cells 
# we need to work a little harder in 2d
nrows = 25 # 13 (50)
ncolumns = 16 # 16
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

im_source_id = sim.add_internal_model(IM_source)
im_id = sim.add_internal_model(IM)

# set all first row cells to have source internal model
sim.set_internal_model(range(ncolumns),im_source_id) # modified
# set all cells but first row to have the same internal model
sim.set_internal_model(range(ncolumns,NCells),im_id) # modified to source also

##########################################################################
# ORDERING BETWEEN ADDING AND SETTING INTERNAL MODELS MUST BE CONSISTENT #
# THIS NEEDS TO BE FIXED                                                 #
##########################################################################

# set up reflecting boundary conditions at bottom edge
# sim.set_boundary_conditions(range(ncolumns),'ref_on')
sim.set_boundary_conditions(range((nrows-1)*ncolumns,nrows*ncolumns),'abs_on')

# need to figure out which cells are connect to which
connections = np.zeros((NCells,NCells)) > 0 # default boolean false array
tri = Delaunay(centers)
for k in xrange(NCells):
    # if index is found in triangle, add indices of points
    # in that triangle to connections
    matches = [t for t in tri.vertices if k in t]
    matches = reduce(lambda x,y: np.concatenate((x,y)),matches)
    # need to filter out connections that are two rows away
    # cells sticking out on the edge get triangulated with cells two rows away
    matches = [m for m in matches if abs(m-k) < 2*ncolumns]
    # duplicates don't really matter since we're just setting them all to True
    connections[k,matches] = True

# use periodic boundary conditions, wrap around
# need to stitch the lattice back together
for i in xrange(nrows):
    # on even rows, the leftmost cell sticks out
    if i % 2 == 0:
        if i > 0:
            # end of last row
            connections[i*ncolumns,i*ncolumns-1] = True
        if i < nrows - 1:
            # end of next row
            connections[i*ncolumns,(i+2)*ncolumns-1] = True
    else:
        if i > 0:
            # end of last row
            connections[(i+1)*ncolumns-1,(i-1)*ncolumns] = True
        if i < nrows - 1:
            # end of next row
            connections[(i+1)*ncolumns-1,(i+1)*ncolumns] = True
    # end of current row
    connections[(i+1)*ncolumns-1,i*ncolumns] = True
    connections[i*ncolumns,(i+1)*ncolumns-1] = True

# for sourcing h, we're going to eliminate diffusion out of the first row
for i in xrange(NCells):
    for j in xrange(NCells):
        if i in range(ncolumns):
            connections[i,j] = False 

# calculate custom square distance matrix for periodic boundary conditions on diffusion
Dmat = np.zeros((NCells,NCells))
for i in xrange(NCells):
    for j in xrange(NCells):
        dx,dy = abs(centers[i,:] - centers[j,:]) # take absolute value to order h to l
        # check if connecting outer boundaries
        lr_connect =  (i % ncolumns == 0 and j % ncolumns == ncolumns - 1)
        rl_connect = (i % ncolumns == ncolumns - 1 and j % ncolumns == 0)
        if lr_connect or rl_connect:
            # we're connecting boundaries
            dx = dx - ncolumns
        Dmat[i,j] = dx ** 2 + dy ** 2

# add diffusion to h and u
# try using custom square distance matrix for period BC
sim.add_interaction('h','h','diffusion',connections,params=([Dh/Th],Dmat))
sim.add_interaction('u','u','diffusion',connections,params=([Du/Tu],Dmat))

# start with only first R8 row up
h_const = 0.015
low_dict = {'a':0.0,'s':0.0,'h':0.0,'u':0.0}
first_low_dict = {'a':0.0,'s':0.0,'h':h_const,'u':0.0} # for first row h
high_dict = {'a':1.0+F,'s':1.0,'h':h_const,'u':0.0} # template row has source h
sim.set_initial_conditions(range(ncolumns,NCells),low_dict)
sim.set_initial_conditions(range(ncolumns),first_low_dict)
sim.set_initial_conditions([4,12],high_dict) # [4,12]

print 'starting simulation'
t = np.linspace(0,250,150)
cdata = sim.simulate(t)
print 'simulation done'


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

# plot 2D picture
def plot_2D(species_name,time):
    tindex = np.abs(t-time).argmin()
    astatus = np.zeros(NCells)
    for i in xrange(NCells):
        astatus[i] = cdata[i][tindex,IM.get_node_id(species_name)]

    # plot using drawn rectangles
    fig = plt.figure()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    # set aspect ratios
    ax1.set(aspect='equal')
    ax2.set(aspect=15)

    # set axis limits
    x,y = centers[:,0],centers[:,1]
    ax1.axis([min(x)-1., max(x)+1., min(y)-1., max(y)+1.])
    
    # data normalizer for color
    data_norm = mpl.colors.Normalize(vmin=min(astatus),vmax=max(astatus))

    dx = [0.8]*NCells # cells are roughly 1 unit apart
    for x,y,c,h in zip(centers[:,0],centers[:,1],plt.cm.jet(data_norm(astatus)),dx):
        ax1.add_artist(Rectangle(xy=(x-h/2.,y-h/2.),facecolor=c,
                                 width=h,height=h,edgecolor='black'))
    
    # now plot colorbar
    mpl.colorbar.ColorbarBase(ax2,cmap=plt.cm.jet,norm=data_norm)

    # plt.scatter(centers[:,0],centers[:,1],
    #             c=astatus,s=100,marker='s',edgecolors='none')
    plt.show()

# we need some better functions for plotting 2D

# plot species at various times
def plot_species(species_name,times,cids=range(NCells)):
    NCells = len(cids)
    tindices = [np.abs(t-v).argmin() for v in times]
    astatus = np.zeros((NCells,len(tindices)))
    for i in xrange(len(cids)):
        for j in xrange(len(tindices)):
            astatus[i,j] = cdata[cids[i]][tindices[j],IM.get_node_id(species_name)]
    '''
    plt.figure()
    plt.scatter(x_coord,np.zeros(NCells),
                c=astatus,s=50,marker='s',edgecolors='none')
    plt.show()
    '''
    # plot along x axis
    plt.figure()
    colors = cm.Dark2(np.linspace(0, 1, len(tindices)))
    x_coord = np.linspace(1,NCells,NCells)
    for j in xrange(len(tindices)):
        plt.scatter(x_coord,astatus[:,j],color=colors[j],s=50)
    plt.legend(['%.1f'% time for time in times],loc='best')
    plt.title(species_name,fontsize=24)
    plt.show()

def plot_all_species(times):
    for n in IM.node_names.keys():
        plot_species(n,times)

def plot_cell(cid,times):
    tstart,tend = [np.abs(t-v).argmin() for v in times]
    plt.figure()
    plt.plot(t[tstart:tend],cdata[cid][tstart:tend,:])
    plt.legend(IM.get_node_names())
    plt.show()

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
    plt.show()

def plot_all_cell_species(time):
    for n in IM.node_names.keys():
        plot_cell_species(n,time)
