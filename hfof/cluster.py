"""
Friends-of-Friends (FOF) for N-body simulations

Peter Creasey - Oct 2016

"""
from __future__ import absolute_import, print_function
from .lib import get_cells, fof3d_periodic, fof_periodic64, get_blocks_cells, minmax
from .primes import smallest_prime_atleast
from numpy import flatnonzero, concatenate, argsort, array, floor, zeros, \
    empty_like, unique, arange
import math


def pad_cube(pos, boxsize, r_pad):
    """
    For a set of points, find images in [0,boxsize)^ndim, then add the repeats
    of those in [0, r_pad) and clone them into [boxsize, boxsize+r_pad)

    An array with the images and the new positions is returned, 
    along with their corresponding indices in the original array (in case you 
    wanted to clone other properties such as their weights).

    pos     - (N,ndim) array
    boxsize - size for periodicity
    r_pad   - edge to pad onto [boxsize, boxsize+r_pad)
    
    returns pad_idx, pos
       pad_idx - (N_new) indices of the positions used to pad
       new_pos - (N+N_new, ndim) array of orig pos + padded pos [0, boxsize+r_pad]
    
    """

    pos = array(pos)
    npts, ndim = pos.shape

    inv_boxsize = float(1.0/boxsize)
    
    scale_r_pad = float(inv_boxsize * r_pad)

    spos = array(pos)*inv_boxsize 
    spos -= floor(spos) # now in [0,1)

    for ax in range(ndim):
        rep_right = zeros((ndim,),spos.dtype)
        rep_right[ax]=1

        # check those within r_pad of 0 
        rt_orig = flatnonzero(spos[:,ax]<scale_r_pad)

        if ax==0:
            # No results from previous padding
            pad_pos = spos[rt_orig]+rep_right
            pad_idx = rt_orig
            continue


        # Some of the *padded* positions may need to be repeated
        rt_pad = flatnonzero(pad_pos[:,ax]<scale_r_pad)

        pad_idx = concatenate((pad_idx, rt_orig, pad_idx[rt_pad]))
            
        pad_pos = concatenate((pad_pos, spos[rt_orig]+rep_right, 
                               pad_pos[rt_pad]+rep_right), axis=0)


    new_pos = concatenate((spos*boxsize, pad_pos*boxsize), axis=0)
                          
    return pad_idx, new_pos

def fof(pos, rcut, boxsize=None, log=None):
    """

    Friends of friends domains (with optional periodicity)

    Return integers for friends-of-friends domains
    """
    assert(pos.shape[1]==3)

    if boxsize is not None:
        old_idx, pos = pad_cube(pos, boxsize, rcut)

    pos_min, pos_max = minmax(pos)
    
    box_dims = pos_max - pos_min
    max_dim = max(box_dims)

    # Match linking length for furthest point in cell
    cell_width = float(rcut / (3**0.5))
    inv_cell_width = float(1.0/cell_width)
    
    if boxsize is not None and boxsize<4*rcut:
        # Cant split into 4x4x4 blocks, have to use old method
        n_min = int(math.ceil(max_dim*inv_cell_width))+2 # two extra cells so never wrap
        
        if log is not None:
            print('Finding primes', file=log)
        N = smallest_prime_atleast(n_min) # Make sure a prime
        M = smallest_prime_atleast(N*N)

        if log is not None:
            print('Searched', N-n_min,'+', M-N*N, 'composites', file=log)
            print('N=%9d (prime index conversion factor)'%N, hex(N), file=log)
            print('M=%9d (prime index conversion factor)'%M, hex(M), file=log)
            
            print('Inserted {:,} images'.format(len(old_idx)), file=log)
            
            print('Position minima', pos_min, file=log)
            print('Position maxima', pos_max, file=log)
            
            
            print('rcut', rcut,file=log)
            
            print('Finding cells', file=log)
        cells = get_cells(pos, inv_cell_width, N, M, log)
    
        if log is not None:
            print('Sorting cells', file=log)        
        sort_idx = argsort(cells)
        if log is not None:
            print('3d fof periodic', file=log)
        domains = fof3d_periodic(cells, N, M, n, old_idx, rcut, sort_idx, pos, log=log)
        return domains


    # Use 4x4x4 block method:
    inv_block_width = 0.25 * inv_cell_width

    # 2 extra blocks so never overlap
    n_min = int(math.ceil(max_dim*inv_block_width))+2

    if boxsize is not None:
        n_min += 1 # 1 block extra for images

    if log is not None:
        print('Finding primes', file=log)
    N = smallest_prime_atleast(n_min) # Make sure a prime
    M = smallest_prime_atleast(N*N)

    if log is not None:
        print('Searched', N-n_min,'+', M-N*N, 'composites', file=log)
        print('N=%9d (prime index conversion factor)'%N, hex(N), file=log)
        print('M=%9d (prime index conversion factor)'%M, hex(M), file=log)

        if boxsize is not None:
            print('Inserted {:,} images'.format(len(old_idx)), file=log)
        
        print('Position minima', pos_min, file=log)
        print('Position maxima', pos_max, file=log)
        
    
        print('rcut', rcut,file=log)
    
        print('Finding cells', file=log)
    blocks_cells = get_blocks_cells(pos, inv_cell_width, N, M, pos_min, log)
    
    if log is not None:
        print('Sorting cells', file=log)        
    sort_idx = argsort(blocks_cells)
    if log is not None:
        print('3d fof', file=log)
    if boxsize is None:
        domains = fof_periodic64(blocks_cells, N, M, rcut, sort_idx, pos, log=log)
    else:
        domains = fof_periodic64(blocks_cells, N, M, rcut, sort_idx, pos, periodic_pad_idx=old_idx, log=log)

    return domains


