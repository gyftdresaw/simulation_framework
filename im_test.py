from Simulation import *

# some testing
IM = InternalModel()
IM.add_node('a','linear',[5])
IM.add_node('b','parabolic',[0.5])
eid = IM.add_edge('a','a','const_prod',params=[2])
IM.add_edge('b',eid,'lin_activ',is_mod=True,mod_type='mult',params=[10])

IM2 = InternalModel()
IM2.add_node('c','parabolic',[2])
IM2.add_edge('c','c','hill_activ',params=[3,1,2])

cell1 = Cell()
cell2 = Cell()
cell3 = Cell()
sim = Simulation()
sim.add_cell(cell1)
sim.add_cell(cell2)
sim.add_cell(cell3)
im_id = sim.add_internal_model(IM)
im2_id = sim.add_internal_model(IM2)

sim.set_internal_model([0,1],im_id)
sim.set_initial_conditions([0],{'a':5,'b':5})
sim.set_initial_conditions([1],{'a':10,'b':10})
sim.set_internal_model([2],im2_id)
sim.set_initial_conditions([2],{'c':0.2})


t = np.linspace(0,10,100)
cdata = sim.simulate(t)
