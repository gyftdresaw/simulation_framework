
# holding external interactions 
import numpy as np
import scipy.spatial.distance as dist

# interaction classes
class Diffusion:
    # connections is a (total_cells x total_cells) boolean matrix 
    # defining which cells interact
    # True for interaction
    def __init__(self,connections,cells,BCs,params=None):
        self.external = True
        self.num_params = 1 # includes constant prefactor
        self.params_set = False
        self.is_mod = False
        self.mod_type = None
        if not params is None:
            self.set_params(params)
        # precalculate distance prefactors
        self.dist_pre = None # (total_cells x total_cells)
        self.set_dist_pre(cells)
        self.cmask = ~connections # in mask, True removes the entry
        self.set_BC_pre(BCs,connections)
    def set_dist_pre(self,cells):
        if self.Dmat_set:
            distances = self.Dmat
        else:
            # should maybe check that all cells have a position set
            # we're going to use scipy distance helpers
            # imported as dist
            pos = np.array([c.position for c in cells])

            # square euclidean distance
            cond_dist = dist.pdist(pos,'sqeuclidean')
            distances = dist.squareform(cond_dist)

        # mask so inversion doesn't break anything
        dmask = (distances == 0.0) # so we don't get that annoying divide by zero
        masked_dists = np.ma.array(distances,mask = dmask,fill_value=0)
        invert_dists = np.ma.divide(1.0,masked_dists)
        self.dist_pre = invert_dists.filled()
        return

    # to be called after set_dist_pre
    # uses BCs to modify dist_pre to mimic 
    # the appropriate boundary conditions
    def set_BC_pre(self,BCs,connections):
        multiplier = np.ones(np.shape(self.dist_pre)) # for reflecting
        self.add_level_diff = np.zeros(np.shape(connections)) # for absorbing
        # add_level_diff indicates which contributions need to be modified for proper
        # absorbing boundary conditions
        for (t,cid_list) in BCs:
            if t == 'ref_on':
                # "reflect on" boundary cells
                # don't multiply other cells on boundary
                dont_mult = np.ones(np.shape(self.dist_pre)) > 0
                for cid in cid_list:
                    dont_mult[cid,:] = False
                    dont_mult[cid,cid_list] = True
                mult_factor = np.ones(np.shape(self.dist_pre))
                mult_factor[(connections * (~dont_mult) > 0)] = 2.0
                multiplier *= mult_factor
            elif t == 'abs_on':
                # absorptive boundary condition
                # "reflect" boundary cells with zero values
                dont_ref = np.ones(np.shape(connections)) > 0
                for cid in cid_list:
                    dont_ref[cid,:] = False
                    dont_ref[cid,cid_list] = True
                self.add_level_diff[(connections * (~dont_ref) > 0)] = 1.0
        self.dist_pre *= multiplier # now dist_pre accounts for BC duplicates
        return

    # normally diffusion will only have one parameter, the diffusion constant
    # to help accomodate periodic boundary conditions or unusual geometry
    # we'll allow for a distance matrix to separately input
    #
    # params is either [Dval] or ([Dval],Dmat)
    def set_params(self,params):
        if len(params) == 1:
            self.C = params[0]
            self.Dmat_set = False
        elif len(params) == 2:
            Dparams = params[0]
            self.Dmat = params[1]
            self.C = Dparams[0]
            self.Dmat_set = True
        self.params_set = True
    # x is an (im_num_cells x total_cells) input matrix
    # im_bounds is a tuple which says which cells this IM corresponds to
    def apply(self,x,im_bounds):
        lower,upper = im_bounds
        # need to extract levels inside cells
        inner_levels = np.diag(x[:,lower:upper])

        # this is kind of a hack
        # should subtract each column elementwise by inner_levels
        level_diff = np.subtract(np.transpose(x),inner_levels).transpose()

        # modify level_diff according to add_level_diff (for absorbing boundary)
        to_add = np.tile(inner_levels,(x.shape[1],1))
        to_add = (-to_add.transpose() * self.add_level_diff[lower:upper,:])

        # multiply through by distance prefactors and constant
        diff_contribs = self.C * self.dist_pre[lower:upper,:] * (level_diff + to_add)
        # apply connection mask
        cmask_slice = self.cmask[lower:upper,:]
        masked_contribs = np.ma.array(diff_contribs,mask=cmask_slice,fill_value=0)

        # finally just sum along the rows
        # filling in 0 for masked entries
        return np.sum(masked_contribs.filled(),axis=1)

# cell-cell activation
class HillActivation:
    # connections is a (total_cells x total_cells) boolean matrix 
    # defining which cells interact
    # True for interaction
    # for now doesn't use cell properties
    def __init__(self,connections,BCs,is_mod=False,mod_type=None,params=None):
        self.external = True
        self.num_params = 1 # includes constant prefactor
        self.params_set = False
        self.is_mod = is_mod
        self.mod_type = mod_type
        if not params is None:
            self.set_params(params)
        self.cmask = ~connections # in mask, True removes the entry
        self.BC_mult = None # (total_cells x total_cells) BC mult factor
        self.set_BC_pre(self,BCs,connections)
    # uses BCs to modify BC_mult to mimic 
    # the appropriate boundary conditions
    def set_BC_pre(self,BCs,connections):
        multiplier = np.ones(np.shape(self.dist_pre))
        for (type,cid_list) in BCs:
            if type == 'ref_on':
                # "reflect on" boundary cells
                # don't multiply other cells on boundary
                dont_mult = np.ones(np.shape(self.dist_pre)) > 0
                for cid in cid_list:
                    dont_mult[cid,:] = False
                    dont_mult[cid,cid_list] = True
                mult_factor = np.ones(np.shape(self.dist_pre))
                mult_factor[(connections * (~dont_mult) > 0)] = 2.0
                multiplier *= mult_factor
        self.BC_mult = multiplier # now dist_pre accounts for BC duplicates
        return
    # C (x/A)^n / (1 + (x/A)^n)
    # params order: [C A n]
    def set_params(self,params):
        self.C = params[0]
        self.A = params[1]
        self.n = params[2]
        self.params_set = True
    # x is an (im_num_cells x total_cells) input matrix
    # im_bounds is a tuple which says which cells this IM corresponds to
    def apply(self,x,im_bounds):
        lower,upper = im_bounds

        # calculate contributions
        hill_contribs = ( self.C * (x/self.A)**self.n / (1 + (x/self.A)**self.n) ) * self.BC_mult[lower:upper,:]
        # apply connection mask
        cmask_slice = self.cmask[lower:upper,:]
        masked_contribs = np.ma.array(hill_contribs,mask=cmask_slice,fill_value=0)

        # finally just sum along the rows
        # filling in 0 for masked entries
        return np.sum(masked_contribs.filled(),axis=1)


# deal with getting the right interaction model
def get_int_model(type,connections,cells,BCs,is_mod=False,mod_type=None,params=None):
    if type == 'diffusion':
        return Diffusion(connections,cells,BCs,params)
    elif type == 'hill_activ':
        return HillActivation(connections,BCs,is_mod,mod_type,params)
    else:
        return None
