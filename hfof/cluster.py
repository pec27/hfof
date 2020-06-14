"""
Friends-of-Friends (FOF) for N-body simulations

Peter Creasey - 2016-2020

"""
from __future__ import absolute_import, print_function
from .lib import get_cells, fof3d_periodic, fof_periodic64, fof_periodic64_2d, \
    get_blocks_cells, get_blocks_cells_2d, minmax, pad_scaled_cube, pad_scaled_square
from .primes import smallest_prime_atleast
from numpy import argsort, array, zeros, arange, float64
import math

def fof(pos, rcut, boxsize=None, log=None):
    """
    Friends of friends domains (with optional periodicity)

    pos  - (N,3) array of positions. For performance this is ideally a numpy C-contiguous array
    rcut - linking length

    [boxsize=None] - periodic box width
    [log=None]     - optional logging

    Return integers for friends-of-friends domains
    """
    assert(pos.shape[1]==3)

    if boxsize is not None:
        if log is not None:
            print('Padding positions', file=log)
        old_idx, pos = pad_scaled_cube(pos, boxsize, rcut, log)

        rcut = rcut / boxsize # scaling applied
        max_dim = 1.0 + rcut # Scaled box + pad
        pos_min = zeros(3, dtype=float64)
    else:
        if log is not None:
            print('Finding minima & maxima', file=log)

        pos_min, pos_max = minmax(pos)
        
        box_dims = pos_max - pos_min
        max_dim = max(box_dims)

    # Match linking length for furthest point in cell
    cell_width = float(rcut / (3**0.5))
    inv_cell_width = float(1.0/cell_width)
    
    if boxsize is not None and rcut>0.25:
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
        else:
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



def fof2d(pos, rcut, boxsize=None, log=None):
    """
    As for fof(..) but a 2-d array of pos

    """
    assert(pos.shape[1]==2)

    if boxsize is not None:
        if log is not None:
            print('Padding positions', file=log)
        old_idx, pos = pad_scaled_square(pos, boxsize, rcut, log)

        rcut = rcut / boxsize # scaling applied
        max_dim = 1.0 + rcut # Scaled box + pad
        pos_min = zeros(2, dtype=float64)
    else:
        if log is not None:
            print('Finding minima & maxima', file=log)

        pos_min, pos_max = minmax(pos)
        
        box_dims = pos_max - pos_min
        max_dim = max(box_dims)

    # Match linking length for furthest point in cell
    cell_width = float(rcut / (2**0.5))
    inv_cell_width = float(1.0/cell_width)
    
    if boxsize is not None and rcut>0.125:
        # Cant split into blocks of 8x8 cells with periodic images because we cant make 1 full block
        
        raise Exception('Periodic box with linking length too large for this method. Usually this means you either \
have very few points (in which case you may want to brute-force), or all points are in a connected group.')



    # Blocks are 8x8 cells
    inv_block_width = 0.125 * inv_cell_width

    # 2 extra blocks so never overlap
    n_min = int(math.ceil(max_dim*inv_block_width))+2

    if boxsize is not None:
        n_min += 1 # 1 block extra for images

    if log is not None:
        print('Finding primes', file=log)
    N = smallest_prime_atleast(n_min) # Make sure a prime

    if log is not None:
        print('Searched', N-n_min, 'composites', file=log)
        print('N=%9d (prime index conversion factor)'%N, hex(N), file=log)

        if boxsize is not None:
            print('Inserted {:,} images'.format(len(old_idx)), file=log)
        else:
            print('Position minima', pos_min, file=log)
            print('Position maxima', pos_max, file=log)
        
    
        print('rcut', rcut,file=log)
    
        print('Finding cells', file=log)

    blocks_cells = get_blocks_cells_2d(pos, inv_cell_width, N, pos_min, log)
    
    if log is not None:
        print('Sorting cells', file=log)        
    sort_idx = argsort(blocks_cells)
    if log is not None:
        print('2d fof', file=log)
    if boxsize is None:
        domains = fof_periodic64_2d(blocks_cells, N, rcut, sort_idx, pos, log=log)
    else:
        domains = fof_periodic64_2d(blocks_cells, N, rcut, sort_idx, pos, periodic_pad_idx=old_idx, log=log)

    return domains



