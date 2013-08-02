
# test modulators (time dependent modulations)
from Simulation import InternalModel,Cell,Simulation
import numpy as np
import matplotlib.pyplot as plt

IM = InternalModel()
IM.add_node('a','linear',[1.0])

a_prod = IM.add_edge('a','a','const_prod',params=[1.0])

cell1 = Cell([0.0])
cell2 = Cell([1.0])
cell3 = Cell([2.0])

sim = Simulation()

cell1_id = sim.add_cell(cell1)
cell2_id = sim.add_cell(cell2)
cell3_id = sim.add_cell(cell3)
IM_id = sim.add_internal_model(IM)

sim.set_internal_model([cell1_id,cell2_id,cell3_id],IM_id)

def lin(x,t):
    return np.maximum((-x+t+np.sin(10.0*t))*(t < 6.0),1.0)

sim.add_modulator(a_prod,'xtgeneral',IM_id,fun=lin)

sim.set_initial_conditions([cell1_id,cell2_id,cell3_id],{'a':0})

t = np.linspace(0,10,1000)
cell_data = sim.simulate(t)

plt.plot(t,cell_data[cell1_id])
plt.plot(t,cell_data[cell2_id])
plt.plot(t,cell_data[cell3_id])
plt.show()
