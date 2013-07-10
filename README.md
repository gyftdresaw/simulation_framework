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

The string `edge_type` specifies the type of interaction and `params` sets the parameters associated with these interactions. Currently the choices of interactions are *'const_prod', 'lin_activ', 'hill_activ', 'hill_inactiv'*; adding new types of interactions is easy and discussed further below. Adding a *'hill_activ'* edge with `params=[C A n]`, corresponds to adding a term *C ([from_node]/A)^n / ( 1 + ([from_node]/A)^n)* to *d[to]/dt*. 

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










