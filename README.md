simulation_framework
====================

This is a Python framework for doing gene network type simulations. After the user specifies what species and species-species interactions exist in what cells, this framework will sort out evaluating and simulating the corresponding differential equations using numpy/scipy.

dependencies
------------



walk through
============

In the main file Simulation.py, there are three central objects: **InternalModel**, **Cell**, and **Simulation**:
 - The `InternalModel` class represents a gene network internal to a cell or a group of cells. Accordingly an InternalModel keeps track of the names of all species in that network as well as all **internal interactions** between species. By 'internal interactions,' I mean all interactions that occur within a single cell. Interactions between species of different cells is handled by the `Simulation` class. 
 - The `Cell` class represents a cell in a simulation. For purposes of simulation, each cell needs an `InternalModel` that describes its inner workings, a set of initial conditions, and any other relevant cell-specific properties defined (position, etc.)
 - The `Simulation` class represents an all-encompassing object that handles evaluating and solving the differential equations associated with a system of cells and their corresponding gene networks. Since a `Simulation` object keeps track of all cells and their internal gene networks, it is natural that **external interactions** between species of different cells are delineated here. 

* * *

We'll now go through a simple 3 species oscillator example (code can be found in negative_feedback.py).

We start by creating a new `InternalModel` and adding our three species to it:

```Python
from Simulation import InternalModel,Cell,Simulation

IM = InternalModel()

IM.add_node('a','linear',[1.0]) # species 'a' decays linearly with constant factor 1.0
IM.add_node('b','linear',[1.0]) # species 'b' does as well
IM.add_node('c','linear',[10.0]) # species 'c' decays 10 times as fast as 'a' or 'b'
```

To add a species (node) to an `InternalModel`, we call `IM.add_node(species_name,degradation=None,params=None)`. Here, `species_name` is some string identifier of the species, `degradation` is a string that specifies the type of degradation is associated with the new species. Currently, the choice is between `'linear'` (goes like C[species]) and `'parabolic'` (goes like C[species]^2). `params` is a list of parameter values associated with the degradation; since only the constant factor *C* needs to be specified, we provide a list with a single value. If `degradation` is set to None, then our new species will not degrade: `IM.add_node('x',None,None)` or `IM.add_node('x')` both add a non-degrading species *x*.

Next, we specify the interactions (edges) amongst the species in our network:

```Python
IM.add_edge('a','a','const_prod',params=[5.0])       # 'a' is produced at a constant rate
IM.add_edge('a','b','hill_activ',params=[5.0,1.0,2]) # 'a' activates 'b'
IM.add_edge('b','c','hill_activ',params=[6.0,1.0,8]) # 'b' activates 'c' (with a sharper cutoff)
```

To add an internal interaction (edge) to an `InternalModel`, we call `IM.add_edge(from_node,to,edge_type,is_mod=False,mod_type=None,params=None)`. All internal interactions originate from some node `from_node` and affect either another node or another edge specifie by `to`. Here, we specify 3 node-node interactions, so both `from_node` and `to` are species string identifiers. 

The string `edge_type` specifies the type of interaction and `params` sets the parameters associated with these interactions. Currently the choices of interactions are *'const_prod', 'lin_activ', 'hill_activ', 'hill_inactiv'*; adding new types of interactions is easy and discussed further below. Adding a *'hill_activ'* edge with `params=[C,A,n]`, corresponds to adding a term *C ([from_node]/A)^n / ( 1 + ([from_node]/A)^n)* to *d[to]/dt*. 

If we want to specify an internal interaction between a node and an edge, we'll need to set the `is_mod` flag to `True` since this interaction is 'modifying' an existing interaction. Before we can add a node->edge interaction, we need some way of identifying the edge being affected. The `InternalModel` class does this by simply assigning an integer edge_id to each edge when added. 

We want to add an interaction in which 'c' inhibits 'a' production. First we'll grab `edge_id` (returned by `IM.add_edge(...)` associated with 'a' production by modifying the last block of code:

```Python
a_prod = IM.add_edge('a','a','const_prod',params=[5.0]) # 'a' is produced at a constant rate
IM.add_edge('a','b','hill_activ',params=[5.0,1.0,2])    # 'a' activates 'b'
IM.add_edge('b','c','hill_activ',params=[6.0,1.0,8])    # 'b' activates 'c' (with a sharper cutoff)
```

Now we can specify an edge ending on `a_prod`:

```Python
IM.add_edge('c',a_prod,'hill_inactiv',is_mod=True,mod_type='mult',params=[1.0,1.0,0.3,8])
```

The flag `mod_type` specifies how this new edge modifies the edge `a_prod`. By setting `mod_type='mult'`, we're saying we want this edge to multiplicatively modify `a_prod`. Specifically, the contribution of `a_prod` to *d['a']/dt* is multiplied by the contribution of this new edge. We could alternatively set `mod_type='intern'` (an *'internal'* modification). In this case, the contribution of this new edge is multiplied into the **input** to the contribution of `a_prod` to *d['a']/dt*. 

In this particular case, the `a_prod` edge contribution to *d['a']/dt* is simply some constant *A*. The contribution of a `'hill_inactiv'` edge with `params=[D,C,A,n]` is *D - C ([from_node]/A)^n / ( 1 + ([from_node]/A)^n)*. Thus, the net contribution to *d['a']/dt* is *A(D - C (['c']/A)^n / ( 1 + (['c']/A)^n))*.


Now we have an internal network specified which schematically looks like *a -> b -> c -| a*. All that's left to do is to make some cells, set some initial conditions and simulate using a `Simulation` object.

Here we make a new `Cell` and `Simulation` object, and then add the `Cell` and `InternalModel` we previously created to the `Simulation` object:

```Python
cell = Cell()                      # new cell with no position specified
sim = Simulation()                 # new Simulation object to organize everything
cell_id = sim.add_cell(cell)       # add the cell to the simulation framework
im_id = sim.add_internal_model(IM) # add the internal model to the simulation framework
```

If we wanted to make a cell with a position (perhaps for diffusion), we would have simply needed to fill in the `position` argument in the call `Cell(position=None)`. We need to add each cell and each internal model we want to simulate to the `Simulation` object; the calls `sim.add_cell(...)` and `sim.add_internal_model(...)` both return corresponding integer id's (starting at 0) for identification. 

After adding the cells and internal models to the `Simulation` object, we need to specify which cells are governed by which internal models and also set the cells initial conditions:

```Python
sim.set_internal_model([cell_id],im_id)
sim.set_initial_conditions([cell_id],{'a':0.1,'b':0.0,'c':0.0})
```

The call `sim.set_internal_model(cell_id_list,im_id)` tells the `Simulation` object that each cell corresponding to the cell id's in `cell_id_list` have the internal model corresponding to `im_id`. Similarly, the call `sim.set_initial_conditions(cell_id_list,ic_dict)` tells the `Simulation` object to set the initial conditions of each cell referenced in `cell_id_list` to those given in the initial condition dictionary `ic_dict`. 

Finally, simulating the system is a simple call:

```Python
import numpy as np

t = np.linspace(0,20,1000)
cell_data = sim.simulate(t)
```

`cell_data` is a list of species time evolution data for each cell. Specifically `cell_data[cell_id]` is a numpy array of dimension (len(t) x num_species) in which rows correspond to the state of each species in that cell's internal model at each time in `t` (1000 time points from 0.0 to 20.0). Plotting this time evolution data is easy:

```Python
import matplotlib.pyplot as plt

plt.plot(t,cell_data[cell_id])
plt.legend(['a','b','c'])
plt.show()
```

So far, we've just been simulating one cell. It is easy to simulate multiple cells with the same internal model and even multiple cells with different internal models:

```Python

# IM1 is the original oscillatory network
IM1 = InternalModel()
# first add the species in the network
IM1.add_node('a','linear',[1.0]) # species 'a' decays linearly with constant factor 1.0
IM1.add_node('b','linear',[1.0]) # species 'b' does as well
IM1.add_node('c','linear',[10.0]) # species 'c' decays 10 times as fast as 'a' or 'b'


a_prod = IM1.add_edge('a','a','const_prod',params=[5.0]) # 'a' is produced at a constant rate
IM1.add_edge('a','b','hill_activ',params=[5.0,1.0,2]) # 'a' activates 'b'
IM1.add_edge('b','c','hill_activ',params=[6.0,1.0,8]) # 'b' activates 'c' (with sharper cutoff)

IM1.add_edge('c',a_prod,'hill_inactiv',is_mod=True,mod_type='mult',params=[1.0,1.0,0.3,8])

# IM2 is a network with only a single species
IM2 = InternalModel()
IM2.add_node('a','linear',[1.0])

# we have 3 cells
cell1 = Cell([0.0])
cell2 = Cell([1.0])
cell3 = Cell([2.0])

sim = Simulation()
# add the cells and internal models to simulation
cell1_id = sim.add_cell(cell1)
cell2_id = sim.add_cell(cell2)
cell3_id = sim.add_cell(cell3)
im1_id = sim.add_internal_model(IM1)
im2_id = sim.add_internal_model(IM2)
```

Here, we create a second internal model that contains the single species 'a' with linear degradation. We also create 3 cells with single dimensional positions 0.0, 1.0, and 2.0, respectively. 

Earlier it was mentioned that *external* interactions (between species of different cells) need to be specified with the `Simulation` object. As an example, we'll allow species 'a' to diffuse between cell2 and cell3:

```Python
# boolean array determines which cells interact
connections = np.array([[False,False,False],[False,True,True],[False,True,True]])
sim.add_interaction('a','a','diffusion',connections,params=[1.0])
```

Analogous to how *internal* interactions were added in an `InternalModel` object, *external* interactions are added in a `Simulation` object with a call to `sim.add_interaction(from_node,to,type,connections,IM_id=None,is_mod=False,mod_type=None,params=None)`. The arguments and flags `from_node`, `to`, `type`, `is_mod`, `mod_type`, and `params` are exactly analogous to those from `IM.add_edge(...)`. The only additional arguments are:
 - connections - (num\_cells x num\_cells) boolean array specifying which cells 
 - IM_id



