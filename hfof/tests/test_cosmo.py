"""
Test data from cosmological simulation (subsampled to 32^3)
"""
from __future__ import absolute_import, print_function
import numpy as np
from hfof import fof, example_data
from numpy.random import RandomState

# linking length:number of domains in the isolated case
ndom_isol = {0.1:2, 0.05:461, 0.01:9600, 0.005:14497, 0.001:27824}

# linking length:number of domains in the periodic case
ndom_prd = {0.1:1, 0.05:396, 0.01:9580, 0.005:14488, 0.001:27823}

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
