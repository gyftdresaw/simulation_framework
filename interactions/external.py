
# holding external interactions 
import numpy as np
import scipy.spatial.distance as dist

# interaction classes
class Diffusion:
    # connections is a (total_cells x total_cells) boolean matrix 
    # defining which cells interact
    # True for interaction
    def __init__(self,connections,cells,params=None):
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
    def set_dist_pre(self,cells):
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
    def set_params(self,params):
        self.C = params[0]
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

        # multiply through by distance prefactors and constant
        diff_contribs = self.C * self.dist_pre[lower:upper,:] * level_diff
        # apply connection mask
        cmask_slice = self.cmask[lower:upper,:]
        masked_contribs = np.ma.array(diff_contribs,mask=cmask_slice,fill_value=0)

        # finally just sum along the rows
        # filling in 0 for masked entries
        return np.sum(masked_contribs.filled(),axis=1)



# deal with getting the right interaction model
def get_int_model(type,connections,cells,is_mod=False,mod_type=None,params=None):
    if type == 'diffusion':
        return Diffusion(connections,cells,params)
    else:
        return None
