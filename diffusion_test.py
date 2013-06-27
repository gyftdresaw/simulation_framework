from Simulation import *

# diffusion testing
IM = InternalModel()
IM.add_node('a') # no degradation

IM2 = InternalModel()
IM2.add_node('a','linear',params=[0.5]) # no degradation
IM2.add_node('b','linear',params=[1])
eid = IM2.add_edge('b','b','hill_activ',params=[1,1,2])
IM2.add_edge('a',eid,'lin_activ',is_mod=True,mod_type='mult',params=[5])

cell1 = Cell([0])
cell2 = Cell([1])
cell3 = Cell([3])

sim = Simulation()
sim.add_cell(cell1)
sim.add_cell(cell2)
sim.add_cell(cell3)

im_id = sim.add_internal_model(IM)
im_id2 = sim.add_internal_model(IM2)

connections = np.array([[True,True,False],[True,True,True],[False,True,True]])

sim.set_internal_model([0,1],im_id)
sim.set_internal_model([2],im_id2)
sim.add_interaction('a','a','diffusion',connections,params=[1])

sim.set_initial_conditions([0],{'a':0})
sim.set_initial_conditions([1],{'a':6})
sim.set_initial_conditions([2],{'a':0,'b':1})

t = np.linspace(0,10,100)
cdata = sim.simulate(t)
