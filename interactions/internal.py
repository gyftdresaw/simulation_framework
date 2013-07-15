
# hold internal interaction classes

import numpy as np

# need constant prefactors
class ConstantProduction:
    def __init__(self,params=None):
        self.num_params = 1 # includes constant prefactor
        self.params_set = False
        self.is_mod = False
        self.mod_type = None
        if not params is None:
            self.set_params(params)
    def set_params(self,params):
        self.C = params[0]
        self.params_set = True
    def apply(self,x):
        return self.C * np.ones(len(x))

class LinearActivation:
    def __init__(self,is_mod=False,mod_type=None,params=None):
        self.num_params = 1 # includes constant prefactor
        self.params_set = False
        self.is_mod = is_mod
        self.mod_type = mod_type
        if not params is None:
            self.set_params(params)
    def set_params(self,params):
        self.C = params[0]
        self.params_set = True
    def apply(self,x):
        return self.C * x

class HillActivation:
    def __init__(self,is_mod=False,mod_type=None,params=None):
        self.num_params = 3 # includes constant prefactor
        self.params_set = False
        self.is_mod = is_mod
        self.mod_type = mod_type
        if not params is None:
            self.set_params(params)
    # C (x/A)^n / (1 + (x/A)^n)
    # params order: [C A n]
    def set_params(self,params):
        self.C = params[0]
        self.A = params[1]
        self.n = params[2]
        self.params_set = True
    def apply(self,x):
        return ( self.C * (x/self.A)**self.n / (1 + (x/self.A)**self.n) )

class HillInactivation:
    def __init__(self,is_mod=False,mod_type=None,params=None):
        self.num_params = 4 # includes constant prefactor
        self.params_set = False
        self.is_mod = is_mod
        self.mod_type = mod_type
        if not params is None:
            self.set_params(params)
    # D - ( C (x/A)^n / (1 + (x/A)^n) )
    # params order: [D C A n]
    def set_params(self,params):
        self.D = params[0]
        self.C = params[1]
        self.A = params[2]
        self.n = params[3]
        self.params_set = True
    def apply(self,x):
        return (self.D - (self.C * (x/self.A)**self.n / (1 + (x/self.A)**self.n) ))

# degradation classes
class LinearDegradation:
    def __init__(self,params=None):
        self.num_params = 1 # includes constant prefactors
        self.params_set = False
        self.is_mod = False
        self.mod_type = None
        if not params is None:
            self.set_params(params)
    def set_params(self,params):
        self.C = params[0]
        self.params_set = True
    # apply assumes given a numpy array
    def apply(self,x):
        return self.C * (-x)

class ParabolicDegradation:
    def __init__(self,params=None):
        self.num_params = 1 # includes constant prefactors
        self.params_set = False
        self.is_mod = False
        self.mod_type = None
        if not params is None:
            self.set_params(params)
    def set_params(self,params):
        self.C = params[0]
        self.params_set = True
    def apply(self,x):
        return self.C * (-(x**2))


# deal with getting the right edge
def get_edge_model(type,is_mod=False,mod_type=None,params=None):
    # interactions
    if type == 'const_prod':
        return ConstantProduction(params)
    elif type == 'lin_activ':
        return LinearActivation(is_mod,mod_type,params)
    elif type == 'hill_activ':
        return HillActivation(is_mod,mod_type,params)
    elif type == 'hill_inactiv':
        return HillInactivation(is_mod,mod_type,params)
    # degradation counts as an interaction
    elif type == 'linear':
        return LinearDegradation(params)
    elif type == 'parabolic':
        return ParabolicDegradation(params)
    # otherwise return None
    else:
        return None
