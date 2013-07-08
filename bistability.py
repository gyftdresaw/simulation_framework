
from Simulation import *
import matplotlib.pyplot as plt

# exclusive bistability test:
# two nodes inhibiting one another

IM = InternalModel()
IM.add_node('a','linear',[0.2])
IM.add_node('b','linear',[0.2])
IM.add_edge('a','b','hill_inactiv',params=[2.0,2.0,2.0,2])
IM.add_edge('b','a','hill_inactiv',params=[2.0,2.0,2.0,2])

cell = Cell()
sim = Simulation()
sim.add_cell(cell)
im_id = sim.add_internal_model(IM)

sim.set_internal_model([0],im_id)
sim.set_initial_conditions([0],{'a':5.0,'b':4.0})

t = np.linspace(0,100,1000)
cdata = sim.simulate(t)

plt.plot(t,cdata[0])
plt.legend(['a','b'])
plt.show()
