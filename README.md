simulation_framework
====================

This is a Python framework for doing gene network type simulations. After the user specifies what species and species-species interactions exist in what cells, this framework will sort out evaluating and simulating the corresponding differential equations using numpy/scipy.

walk through
============

In the main file Simulation.py, there are three central objects: **InternalModel**, **Cell**, and **Simulation**:
 - The **InternalModel** class represents a gene network internal to a cell or a group of cells. Accordingly an InternalModel keeps track of the names of all species in that network as well as all **internal interactions** between species. By 'internal interactions,' I mean all interactions that occur within a single cell. Interactions between species of different cells is handled by the **Simulation** class. 
 - The **Cell** class represents a cell in a simulation. Each cell has an In
