from Simulation import *
from scipy.spatial import Delaunay
import utils.lattice as lattice
import utils.plotting as plotting
import utils.saving as saving
import os

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
nrows = 5 # 13 (50)
ncolumns = 16 # 16
NCells = nrows * ncolumns

# 'center' of each cell on hexagonal lattice
centers = lattice.get_centers(nrows,ncolumns)

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
connections = lattice.get_connections(nrows,ncolumns,periodic=True)
# for sourcing h, we're going to eliminate diffusion out of the first row
for i in xrange(NCells):
    for j in xrange(NCells):
        if i in range(ncolumns):
            connections[i,j] = False 

# calculate custom square distance matrix for periodic boundary conditions on diffusion
Dmat = lattice.get_periodic_dmat(nrows,ncolumns)

# add diffusion to h and u
# try using custom square distance matrix for period BC
sim.add_interaction('h','h','diffusion',connections,params=([Dh/Th],Dmat))
sim.add_interaction('u','u','diffusion',connections,params=([Du/Tu],Dmat))

# start with only first R8 row up
h_const = 0.016
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

# save to pickle file with same name as script
pickle_file = os.path.basename(__file__).split('.')[0] + '.p'
saving.save(pickle_file,t=t,cdata=cdata,centers=centers,sim=sim)

# reload for plotting
PH = plotting.PlotHelper.from_file(pickle_file)

