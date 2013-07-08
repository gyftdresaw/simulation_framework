from Simulation import *
import matplotlib.pyplot as plt

# negative-feedback loop
# should be able to see oscillations, damped and sustained

IM = InternalModel()
IM.add_node('a','linear',
