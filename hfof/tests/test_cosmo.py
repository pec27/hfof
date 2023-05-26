"""
Test data from cosmological simulation (subsampled to 32^3)
"""
from __future__ import absolute_import, print_function
import numpy as np
from hfof import fof, fof2d, example_data
from numpy.random import RandomState

# linking length:number of domains in the isolated case
ndom_isol = {0.1:2, 0.05:461, 0.01:9600, 0.005:14497, 0.001:27824}

# linking length:number of domains in the periodic case
ndom_prd = {0.3: 1, 0.1:1, 0.05:396, 0.01:9580, 0.005:14488, 0.001:27823}

def test_cosmo_isol():
    """ Cosmological (non-periodic) """
    pos = example_data.get_pos()

    for rcut, n in ndom_isol.items():
        domains = fof(pos, rcut)
        n_doms = len(np.unique(domains))
        assert(n==n_doms)
        print(rcut, n_doms, n)


def test_cosmo_periodic():
    """ Cosmological (periodic) """
    pos = example_data.get_pos()

    for rcut, n in ndom_prd.items():
        domains = fof(pos, rcut, boxsize=1.0)
        n_doms = len(np.unique(domains))
        assert(n==n_doms)
        print(rcut, n_doms, n)
    
              
def test_shift_scale():
    """ Cube scaled and shifted """
    shift, scale = (8.82334, 4.234, -15.234), 4.8271

    pos = example_data.get_pos()*scale + shift

    # scaled the box => scale the linking length

    # Isolated
    for rcut, n in ndom_isol.items():
        domains = fof(pos, rcut*scale)
        n_doms = len(np.unique(domains))
        assert(n==n_doms)
        print(rcut, n_doms, n)

    #Periodic (scale boxsize too)
    for rcut, n in ndom_prd.items():
        domains = fof(pos, rcut*scale, boxsize=scale)
        n_doms = len(np.unique(domains))
        assert(n==n_doms)
        print(rcut, n_doms, n)

def test_consistent_groups():
    """ Shuffling gives consistent groups """
    # Note we dont actually need to check for renumbering between hfof and 
    # itself because even if you shuffle the points they end up in the same 
    # cells, which are treated in order, so the group numbers (at each 
    # position) are guaranteed to be the same. Only if you check against 
    # another FOF do you need to check the numbering.

    pos = example_data.get_pos()
    N = pos.shape[0]
    rs = RandomState(seed=123)
    idx = np.arange(N)
    rs.shuffle(idx)
    rcut = 0.01


    # Non=periodic
    doms1 = fof(pos, rcut)
    doms2 = fof(pos[idx], rcut)

    assert(np.all(np.equal(doms1[idx], doms2)))
    # Periodic
    doms1 = fof(pos, rcut, boxsize=0.8)
    doms2 = fof(pos[idx], rcut, boxsize=0.8)

    assert(np.all(np.equal(doms1[idx], doms2)))

def test_cosmo_float32():
    """ Cosmological (non-periodic) """
    pos = np.float32(example_data.get_pos())

    for rcut, n in ndom_isol.items():
        domains = fof(pos, rcut)
        n_doms = len(np.unique(domains))
        assert(n==n_doms)
        print(rcut, n_doms, n)

def test_2d_domains():
    """ 2d domains equivalent to 3d code with 0 in z-dimension """
    pos = example_data.get_pos() # (32768,3) positions in [0,1]
    pos[:,2] = 0.0
    pos2d = np.require(pos[:,:2], requirements=['C'])
    
    for boxsize in [None, 1.0]:
        for r_cut in [0.001, 0.005, 0.1]: # linking length    

            dom3d = fof_3d_labels = fof(pos, r_cut, boxsize=boxsize) # integer labels from 0...

            dom2d = fof_2d_labels = fof2d(pos2d, r_cut, boxsize=boxsize)
            num_dom3d = len(np.unique(dom3d))
            num_dom2d = len(np.unique(dom2d))
            print('Number of domains 3d: {:,}, 2d: {:,}'.format(num_dom3d, num_dom2d))
            assert(num_dom3d==num_dom2d)
