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
IM.add_edge('b','c','hill_activ',params=[6.0,1.0,8]) # 'b' activates 'c'
```







