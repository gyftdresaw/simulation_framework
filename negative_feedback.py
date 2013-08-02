from Simulation import InternalModel,Cell,Simulation
import matplotlib.pyplot as plt
import numpy as np

# negative-feedback loop
# should be able to see oscillations, damped and sustained
# need an intermediate species ('b') in order to see sustained oscillations
# -- the necessary delay can't be achieved with just 'a' and 'c'

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

sim.set_internal_model([cell1_id,cell2_id],im1_id) # first two cells have IM1
sim.set_internal_model([cell3_id],im2_id)          # last cell as IM2

connections = np.array([[False,False,False],[False,True,True],[False,True,True]])
sim.add_interaction('a','a','diffusion',connections,params=[1.0])

# cell1 and cell2 start with the same initial conditions
sim.set_initial_conditions([cell1_id,cell2_id],{'a':0.1,'b':0.0,'c':0.0})
sim.set_initial_conditions([cell3_id],{'a':0}) # cell3 starts with 0 'a'

t = np.linspace(0,20,1000)
cell_data = sim.simulate(t) # cell_data is a list of 3 numpy arrays containing cell specific data

plt.figure()
plt.subplot(311)
plt.plot(t,cell_data[cell1_id])
plt.title('cell1')
plt.legend(['a','b','c'])

plt.subplot(312)
plt.plot(t,cell_data[cell2_id])
plt.title('cell2')
plt.legend(['a','b','c'])

plt.subplot(313)
plt.plot(t,cell_data[cell3_id])
plt.title('cell3')
plt.legend(['a'])
plt.xlabel('time')

plt.tight_layout()
plt.show()
