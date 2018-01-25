"""
Friends-of-Friends (FOF) for N-body simulations

Peter Creasey - Oct 2016

"""
from __future__ import absolute_import, print_function
from hfof.lib import fof3d, get_cells
from lizard.periodic import pad_unitcube
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix, csgraph
from numpy import square, flatnonzero, ones, zeros_like, cumsum, concatenate, \
    arange, searchsorted, bincount, sort, diff, int8, argsort, array, unique, empty, int32, uint32
import sys
import math

def fof(pos, rcut, log=None):
    """
    Return integers for friends-of-friends domains
    """
    n, dim = pos.shape
    assert(dim==3)

    pos_min = pos.min(axis=0)
    pos_max = pos.max(axis=0)
    
    box_dims = pos_max - pos_min
    max_dim = max(box_dims)
    
    
    cell_width = float(rcut / (3**0.5))
    
    inv_cell_width = float(1.0/cell_width)
    
    nval = int(math.ceil(max_dim*inv_cell_width))+2
    nval = nval | 3 # Make sure last two bits set to 1
    
    if log is not None:
        print('Nval=%d (index conversion factor)'%nval, bin(nval), file=log)
        
        print('Position minima', pos_min, file=log)
        print('Position maxima', pos_max, file=log)
        
    
        print('rcut', rcut,file=log)
    
        print('Finding cells', file=log)
    cells = get_cells(pos, inv_cell_width, nval, log)
    
    if log is not None:
        print('Sorting cells', file=log)        
    sort_idx = argsort(cells)
    if log is not None:
        print('3d fof', file=log)
    domains = fof3d(cells, nval, rcut, sort_idx, pos, log=log)
    return domains
