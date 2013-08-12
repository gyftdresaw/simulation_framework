
# Framework for cell-network simulations using scipy/numpy
import numpy as np
from scipy.integrate import odeint
import interactions.internal as intern
import interactions.external as extern
import interactions.modulators as moduls

# need time scale
class Node:
    def __init__(self,name):
        self.name = name
        self.degrades = False # no degradation by default
        self.degradation = None
        self.diffuses = False # since diffusion involves multiple cells
        self.diffusion = None # it should be set with cell-cell interactions
        self.is_primary = False # primary false by default

    def set_degradation(self,degradation):
        self.degrades = not degradation is None
        self.degradation = degradation
        return self.degradation            

    def set_as_primary(self):
        self.is_primary = True

class Edge:
    def __init__(self,from_node,to,model,edge_id):
        self.from_node = from_node
        # edge can either point node -> node
        # or node -> edge (is_mod = True)
        if model.is_mod:
            self.to_node = None
            self.to_edge = to
        else:
            self.to_node = to
            self.to_edge = None
        self.is_mod = model.is_mod
        self.mod_type = model.mod_type # 'intern' or 'mult'
        self.model = model # model contains info relevant to interaction
        self.edge_id = edge_id # simple but not the greatest way to identify edges
        self.external = False

class InternalModel:
    def __init__(self):
        self.nodes = [] # list of species nodes
        self.edges = [] # list of edges
        self.num_params = [] # number of params associated with each node
        self.node_names = {} # dict of name -> node index
        self.edge_counter = 0 # for generating edge ids 

    def get_node(self,name):
        return self.nodes[self.node_names[name]]
    
    def get_node_id(self,name):
        return self.node_names[name]

    def set_node_degradation(self,name,degradation,params=None):
        # get node
        node = self.get_node(name)
        deg_model = intern.get_edge_model(degradation,params=params)
        node.set_degradation(deg_model)

        # also add degradation as an edge
        eid = self.edge_counter
        self.edge_counter += 1
        new_edge = Edge(name,name,deg_model,eid)
        self.edges.append(new_edge)
        return self.get_edge(eid)        

    def add_node(self,name,degradation=None,params=None):
        new_node = Node(name)
        # add node to list to keep track of it
        self.nodes.append(new_node)
        self.node_names[name] = len(self.nodes) - 1 # set name in name dict

        # set degradation of new node
        # call helper
        self.set_node_degradation(name,degradation,params)
        

        return self.get_node(name)

    def get_edge(self,edge_id):
        # simply look for edge_id in edge list
        # take first element of match
        return ([e for e in self.edges if e.edge_id == edge_id])[0]

    def add_edge(self,from_node,to,type,is_mod=False,mod_type=None,params=None):
        eid = self.edge_counter
        self.edge_counter += 1
        new_model = intern.get_edge_model(type,is_mod,mod_type,params)
        new_edge = Edge(from_node,to,new_model,eid)
        self.edges.append(new_edge)
        return eid

    def num_nodes(self):
        return len(self.nodes)

    def get_node_names(self):
        return [n.name for n in self.nodes]

    # helper to make sure we've set all our parameters
    def all_params_set(self):
        # degradation are included in edges
        edges_set = [e.model.params_set for e in self.edges]
        
        # shortcut for and-ing together params_set
        if len(edges_set) > 0:
            return reduce(lambda x,y: x and y, edges_set)
        else:
            return True
    
    # return all edges that contribute to a nodes dynamics
    def get_contributions(self,to_node):
        return [e for e in self.edges if e.to_node == to_node]

    # return all edges that modify this interaction
    def get_modifiers(self,to_edge):
        return [e for e in self.edges if e.to_edge == to_edge and e.is_mod]

    def set_primary_node(self,name):
        return
    
    # return list of node names in id order
    def get_node_names(self):
        return [n.name for n in self.nodes]
    
    # TODO: this would probably be useful
    def print_model(self):
        return 

class Cell:
    def __init__(self,position=None):
        self.position = None
        if not position is None:
            self.set_pos(position)
        self.IM = None # internal model for cell
        self.IC = None # initial conditions 
        self.num_species = 0 # number of dynamic species in cell
        self.cell_id = None

    def set_cell_id(self,cid):
        self.cell_id = cid

    # pos should be vector [x,y] or [x,y,z]
    # might make 1D, 2D or 3D easier later
    def set_pos(self,pos):
        self.position = pos

    def set_IM(self,IM):
        self.IM = IM
        self.num_species = IM.num_nodes()
        return

    # IC is a dictionary
    # but we're going to set the internal boundary conditions
    # according to the ordering given by IM.nodes
    def set_IC(self,IC):
        # check some consistency things
        # make sure we already have an internal model
        if self.IM is None:
            print 'We can\'t set initial conditions without a model'
            return
        if not len(IC) == self.IM.num_nodes():
            print 'We have %d ICs but %d nodes' % (len(IC),self.IM.num_nodes())
            return
        # if we've made it here should be good
        # unless theres an incorrect key
        self.IC = np.zeros(len(IC))
        for i in xrange(len(IC)):
            self.IC[i] = IC[self.IM.nodes[i].name]
        return 

# wrapper for interactions
class Interaction:
    # IM_id only relevant if interaction is a mod of an internal edge
    def __init__(self,from_node,to,model,int_id,IM_id=None):
        self.external = True
        self.from_node = from_node
        # either (node -> node) (is_mod is False)
        # or (node -> edge) -- where edge can be internal or external
        if model.is_mod:
            self.to_node = None
            self.to_edge = to
            self.to_IM = IM_id # if IM_id is None, then refers to interaction
        else:
            self.to_node = to
            self.to_edge = None
            self.to_IM = None

        self.is_mod = model.is_mod
        self.mod_type = model.mod_type # 'intern' or 'mult'

        self.model = model
        self.int_id = int_id

class Modulator:
    # IM_id only relevant if modulator of an internal edge
    def __init__(self,to_edge,model,IM_id=None):
        self.to_edge = to_edge
        self.to_IM = IM_id # if IM_id is None, then modulating an external interaction
        self.model = model
    
class Simulation:
    def __init__(self):
        self.cells = [] # cell list
        self.IMs = [] # list of internal models for cells
        self.IM_cells = [] # list of list of cell ids associated with each IM
        self.IM_bounds = [] # list of tuples denoting cell bounds
        self.interactions = [] # list of cell-cell interactions
        self.int_counter = 0 # for generating int_id's
        self.modulators = [] # list of time-dependent modulators
        
    def add_cell(self,cell):
        self.cells.append(cell)
        cid = len(self.cells) - 1
        cell.set_cell_id(cid)
        return cid # return cell id

    def add_internal_model(self,IM):
        self.IMs.append(IM)
        im_id = len(self.IMs) - 1
        self.IM_cells.append([]) # append empty list of cells associated with IM
        return im_id # return IM id

    def set_internal_model(self,cell_id_list,IM_id):
        for c in [self.cells[cid] for cid in cell_id_list]:
            c.set_IM(self.IMs[IM_id])
        # this is pretty dangerous, should really be checking
        # for multiple assignments, duplicates, etc.
        self.IM_cells[IM_id].extend(cell_id_list) # just append all cell ids to IM list
        return 

    def set_initial_conditions(self,cell_id_list,IC):
        # IC is species dictionary
        for c in [self.cells[cid] for cid in cell_id_list]:
            c.set_IC(IC)
        return

    # need to set constants for each 
    # it would probably be easier to do this with the internal model
    # but might be better (for redundancy and simplification) to do here
    # will do with internal nodes for now
    # def set_prefactors(self,IM):
    #     return
              
    def add_interaction(self,from_node,to,type,connections,IM_id=None,is_mod=False,mod_type=None,params=None):
        int_id = self.int_counter
        self.int_counter += 1
        new_model = extern.get_int_model(type,connections,self.cells,is_mod,mod_type,params)
        new_interaction = Interaction(from_node,to,new_model,int_id,IM_id)
        self.interactions.append(new_interaction)
        return int_id

    def add_modulator(self,to,type,IM_id=None,**kwargs):
        new_model = moduls.get_mod_model(type,self.cells,**kwargs)
        new_modulator = Modulator(to,new_model,IM_id)
        self.modulators.append(new_modulator)
        return

    def num_cells(self):
        return len(self.cells)
    
    # checks:
    # all cells have internal model
    # all cell shave initial conditions set
    # all internal models have parameters set
    def simulation_ready(self):
        cIM_ready_list = [not c.IM is None for c in self.cells]
        cells_IM_ready = reduce(lambda x,y: x and y, cIM_ready_list)
        if not cells_IM_ready:
            print 'some internal models are not set!'
            return False
        cIC_ready_list = [not c.IC is None for c in self.cells]
        cells_IC_ready = reduce(lambda x,y: x and y, cIC_ready_list)
        if not cells_IC_ready:
            print 'some initial conditions are not set!'
            return False
        IMp_ready_list = [im.all_params_set() for im in self.IMs]
        IM_params_ready = reduce(lambda x,y: x and y, IMp_ready_list)
        if not IM_params_ready:
            print 'some parameters not set in an internal model!'
            return False
        # otherwise we're good to go
        return True

    # set IM_bounds
    def set_cell_bounds(self):
        ncells = [len(x) for x in self.IM_cells]
        bounds = []
        low = 0
        for i in xrange(len(ncells)):
            bounds.append((low,ncells[i]+low))
            low += ncells[i]
        self.IM_bounds = bounds
        return


    # return all external interactions that affect this node
    def get_contributions(self,to_node):
        return [i for i in self.interactions if i.to_node == to_node]
    
    # if IM_id is given, then we're looking for modifiers to an internal edge
    # if not, then we're looking for modifiers to an external edge
    def get_modifiers(self,to_edge,IM_id = None):
        return [i for i in self.interactions if i.to_edge==to_edge and i.to_IM==IM_id]

    # if IM_id is given, then we're looking for modulators to an internal edge
    # if not, then we're looking fro modulators to an external edge
    def get_modulators(self,to_edge,IM_id = None):
        return [m for m in self.modulators if m.to_edge==to_edge and m.to_IM==IM_id]

    # recursively resolve the contribution 
    # IM_id is relevant, tells us for which cells we're calculating for 
    def resolve_contribution(self,contrib,ycurrent,IM_id,t):

        # contrib is either an internal or external edge
        # applying them will require different data depending
        if contrib.external:
            # only external modifiers possible on external contributions
            ext_mods = self.get_modifiers(contrib.int_id) # no IM_Id needed
            
            # need to construct xcontrib
            # which will be (IM_num_cells x total_num_cells)
            # note: np.append copies the array
            xdata = np.array([])
            for i in xrange(len(self.IMs)):
                if contrib.from_node in self.IMs[i].node_names:
                    from_id = self.IMs[i].get_node_id(contrib.from_node)
                    xdata = np.append(xdata,ycurrent[i][:,from_id])
                else:
                    # if IM doesn't have that species fill with nans
                    nan_fill = np.empty(len(self.IM_cells[i]))
                    nan_fill.fill(np.nan)
                    xdata = np.append(xdata,nan_fill)

            IM_num_cells = len(self.IM_cells[IM_id]) 
            total_cells = self.num_cells()
            xcontrib = np.tile(xdata,(IM_num_cells,1))

            # mods to xcontrib
            for mod in ext_mods:
                if mod.mod_type == 'intern':
                    # multiply rows of xcontrib through by mod_contrib
                    mod_contrib = resolve_contribution(mod,ycurrent,IM_id,t)
                    xcontrib = np.dot(np.diag(mod_contrib),xcontrib)

            # apply xcontrib, then multiply by multiplicative mods
            # for application need proper bounds
            ycontrib = contrib.model.apply(xcontrib,self.IM_bounds[IM_id])
            for mod in ext_mods:
                if mod.mod_type == 'mult':
                    ycontrib *= self.resolve_contribution(mod,ycurrent,IM_id,t)

            # check for any modulators (time dependent mods)
            modulates = self.get_modulators(contrib.edge_id) # no IM_id needed
            for modulate in modulates:
                ycontrib *= modulate.model.apply(self.IM_cells[IM_id],t)
            return ycontrib
        else:
            # need to check for internal and external modifiers
            int_mods = self.IMs[IM_id].get_modifiers(contrib.edge_id)
            ext_mods = self.get_modifiers(contrib.edge_id,IM_id)
            
            from_node = self.IMs[IM_id].get_node_id(contrib.from_node)
            
            # THIS SLICING DIDN'T MAKE A COPY?!?!
            # arrays generated by basic slicing are VIEWS of the original array
            # unbelievable
            xcontrib = np.copy(ycurrent[IM_id][:,from_node])
            
            # first do mods to xcontrib
            for mod in (int_mods + ext_mods):
                if mod.mod_type == 'intern':
                    xcontrib *= self.resolve_contribution(mod,ycurrent,IM_id,t)
                    
            # apply xcontrib, then multiply by multiplicative mods
            ycontrib = contrib.model.apply(xcontrib)
            for mod in (int_mods + ext_mods):
                if mod.mod_type == 'mult':
                    ycontrib *= self.resolve_contribution(mod,ycurrent,IM_id,t)

            # print from_node,IM_id,contrib.edge_id,ycontrib,ycurrent
            
            # check for any modulators (time dependent mods)
            modulates = self.get_modulators(contrib.edge_id,IM_id)
            for modulate in modulates:
                ycontrib *= modulate.model.apply(self.IM_cells[IM_id],t)
            return ycontrib


    # derivative for odeint
    # y is flattened state of all species of all cells 
    def deriv(self,y,t):
        # print 'getting deriv'

        # need a reshape for each IM: cells x species
        IM_ncells = [len(self.IM_cells[i]) for i in xrange(len(self.IMs))]
        IM_nspecies = [im.num_nodes() for im in self.IMs]
        IM_dims = zip(IM_ncells,IM_nspecies) # list of tuples (ncells,nspecies)
        
        # for extracting data per IM
        bounds = [0] + [a*b for (a,b) in IM_dims]
        # cumulative sum to make bounds correct
        for i in xrange(1,len(bounds)):
            bounds[i] = bounds[i] + bounds[i-1]
        
        # extract data for each IM
        ycurrent  = []
        yprime = []
        for i in xrange(len(IM_dims)):
            ycurrent.append(y[bounds[i]:bounds[i+1]].reshape(IM_dims[i]))
            yprime.append(np.zeros(IM_dims[i]))

        # ycurrent is a list of numpy arrays corresponding to cell data
        # for each internal model
        # and yprime will be the corresponding data for the derivative
        
        # for each internal model, calculate derivative for each species
        for i in xrange(len(self.IMs)):
            # for identifying nodes in this model
            cnodes = self.IMs[i].get_node_names()

            # for each species in this internal model
            for j in xrange(self.IMs[i].num_nodes()):
                # find all interactions that end on this node
                # there can be two sources: internal and external
                int_contributions = self.IMs[i].get_contributions(cnodes[j])
                ext_contributions = self.get_contributions(cnodes[j]) # need this

                # add up all contributions
                # resolve contribution will deal with modifier business
                for contrib in (int_contributions + ext_contributions):
                    yprime[i][:,j] += self.resolve_contribution(contrib,ycurrent,i,t)
                    # print yprime[i][:,j],i,j,ycurrent[i][:,j]

        # print 'done calculating derivative'
        '''
        for i in xrange(len(self.IMs)):
            for j in xrange(self.IMs[i].num_nodes()):
                print ycurrent[i][:,j],i,j
        print 'this is ycurrent'
        '''
        # now just need to reshape yprime and return it
        flat_yprime = [yp.flatten() for yp in yprime]
        # concat the flattened yprimes for each internal model
        return reduce(lambda x,y: np.append(x,y),flat_yprime)

    def simulate(self,t):
        # make sure checks are satisfied
        if not self.simulation_ready():
            return

        # set cell bounds
        self.set_cell_bounds()
        # need to get initial conditions together, and the ordering right
        yinitial = [] # array of (IM_cell x num_species)
        for i in xrange(len(self.IMs)):
            yinitial.append(np.empty((len(self.IM_cells[i]),self.IMs[i].num_nodes())))
            for j in xrange(len(self.IM_cells[i])):
                yinitial[i][j,:] = self.cells[self.IM_cells[i][j]].IC

        # flatten yinitial
        flat_yinitial = [yi.flatten() for yi in yinitial]
        y0 = reduce(lambda x,y: np.append(x,y),flat_yinitial)

        # call odeint
        sol,info = odeint(self.deriv,y0,t,full_output=True)
        
        # deal with results
        # for each cell, have time data (time_points x num_species)
        cdata = [None for c in self.cells]
        
        IM_ncells = [len(self.IM_cells[i]) for i in xrange(len(self.IMs))]
        IM_nspecies = [im.num_nodes() for im in self.IMs]
        IM_dims = zip(IM_ncells,IM_nspecies) # list of tuples (ncells,nspecies)

        sol_index = 0
        for i in xrange(len(self.IMs)):
            for j in xrange(len(self.IM_cells[i])):
                cid = self.IM_cells[i][j]
                cdata[cid] = np.empty((len(t),IM_nspecies[i]))
                cdata[cid] = sol[:,sol_index:(sol_index+IM_nspecies[i])]
                sol_index += IM_nspecies[i]

        return cdata

