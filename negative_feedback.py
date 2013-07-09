from Simulation import *
import matplotlib.pyplot as plt

# negative-feedback loop
# should be able to see oscillations, damped and sustained
# need an intermediate species ('c') in order to see sustained oscillations
# -- the necessary delay can't be achieved with just 'a' and 'b'

IM = InternalModel()
IM.add_node('a','linear',[0.1])
IM.add_node('b','linear',[0.9])
IM.add_node('c','linear',[0.1])
a_self = IM.add_edge('a','a','const_prod',params=[1.0])
IM.add_edge('a','c','hill_activ',params=[5.0,1.0,2])
IM.add_edge('c','b','hill_activ',params=[1.0,1.0,4])
IM.add_edge('b',a_self,'hill_inactiv',is_mod=True,mod_type='mult',params=[1.0,1.0,0.3,8])

cell = Cell()
sim = Simulation()
sim.add_cell(cell)
im_id = sim.add_internal_model(IM)

sim.set_internal_model([0],im_id)
sim.set_initial_conditions([0],{'a':1.0,'b':0.0,'c':0.0})

t = np.linspace(0,300,1000)
cdata = sim.simulate(t)


plt.plot(t,cdata[0])
plt.legend(['a','b','c'])
plt.show()
