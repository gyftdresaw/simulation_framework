
# modulator classes
import numpy as np

class XTGeneral:
    def __init__(self,cells,fun=None):
        self.xvals = np.array([c.position[0] for c in cells])
        self.fun = fun
    def apply(self,cids,t):
        return self.fun(self.xvals[cids],t)

class ZTGeneral:
    def __init__(self,cells,fun=None):
        self.posvals = np.array([c.position for c in cells])
        self.fun = fun
    def apply(self,cids,t):
        return self.fun(self.posvals,t)

def get_mod_model(type,cells,**kwargs):
    if type == 'xtgeneral':
        return XTGeneral(cells,**kwargs)
    elif type == 'ztgeneral':
        return ZTGeneral(cells,**kwargs)
    else:
        return None
