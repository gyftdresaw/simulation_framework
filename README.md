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

We'll now go through a simple example:

We start by 





