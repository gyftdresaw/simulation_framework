
# some convenience functions for creating hexagonal lattices
import numpy as np
from scipy.spatial import Delaunay

# 'center' of each cell in hexagonal lattice
# our indices i and j go like:
# i: left -> right
# j: bottom -> top
def get_centers(nrows,ncolumns):
    NCells = nrows * ncolumns
    centers = np.empty((NCells,2)) # a center for each cell
    for i in xrange(ncolumns):
        for j in xrange(nrows):
            c = np.array([float(i) + 0.5*(j % 2),(np.sqrt(3)/2)*float(j)])
            centers[i + j*ncolumns,:] = c
    return centers


# get connections for hexagonal lattice
def get_connections(nrows,ncolumns,periodic=False):
    NCells = nrows * ncolumns
    connections = np.zeros((NCells,NCells)) > 0 # default boolean false array

    centers = get_centers(nrows,ncolumns)
    tri = Delaunay(centers)
    for k in xrange(NCells):
        # if index is found in triangle, add indices of points
        # in that triangle to connections
        matches = [t for t in tri.vertices if k in t]
        matches = reduce(lambda x,y: np.concatenate((x,y)),matches)
        # need to filter out connections that are two rows away
        # cells sticking out on the edge get triangulated with cells two rows away
        matches = [m for m in matches if abs(m-k) < 2*ncolumns]
        # duplicates don't really matter since we're just setting them all to True
        connections[k,matches] = True
        
    # if using periodic boundary conditions, wrap around
    # need to stitch the lattice back together
    if periodic:
        for i in xrange(nrows):
            # on even rows, the leftmost cell sticks out
            if i % 2 == 0:
                if i > 0:
                    # end of last row
                    connections[i*ncolumns,i*ncolumns-1] = True
                if i < nrows - 1:
                    # end of next row
                    connections[i*ncolumns,(i+2)*ncolumns-1] = True
            else:
                if i > 0:
                    # end of last row
                    connections[(i+1)*ncolumns-1,(i-1)*ncolumns] = True
                if i < nrows - 1:
                    # end of next row
                    connections[(i+1)*ncolumns-1,(i+1)*ncolumns] = True
            # end of current row
            connections[(i+1)*ncolumns-1,i*ncolumns] = True
            connections[i*ncolumns,(i+1)*ncolumns-1] = True
    
    return connections

# calculate custom square distance matrix for periodic boundary conditions on diffusion
def get_periodic_dmat(nrows,ncolumns):
    NCells = nrows * ncolumns
    centers = get_centers(nrows,ncolumns)
    Dmat = np.zeros((NCells,NCells))
    for i in xrange(NCells):
        for j in xrange(NCells):
            dx,dy = abs(centers[i,:] - centers[j,:]) # take absolute value to order h to l
            # check if connecting outer boundaries
            lr_connect =  (i % ncolumns == 0 and j % ncolumns == ncolumns - 1)
            rl_connect = (i % ncolumns == ncolumns - 1 and j % ncolumns == 0)
            if lr_connect or rl_connect:
                # we're connecting boundaries
                dx = dx - ncolumns
            Dmat[i,j] = dx ** 2 + dy ** 2
    return Dmat
