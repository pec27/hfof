"""
Test data from cosmological simulation (subsampled to 32^3)
"""
from __future__ import absolute_import, print_function
import numpy as np
from hfof import fof
from os import path

_data = None
# linking length:number of domains in the isolated case
ndom_isol = {0.1:2, 0.05:461, 0.01:9600, 0.005:14497, 0.001:27824}

# linking length:number of domains in the periodic case
ndom_prd = {0.1:1, 0.05:396, 0.01:9580, 0.005:14488, 0.001:27823}


def get_cosmo_pos():
    """
    Load cosmo simulation data (unit cube, (32768,3) array)
    """
    global _data
    if _data is not None:
        return _data

    name = path.join(path.dirname(path.abspath(__file__)), 'cosmo32768.dat')
    if not path.exists(name):
        raise Exception('File '+str(name)+' does not exist. Maybe I dont understand your OS?')

    pos = np.reshape(np.fromfile(name, np.uint16, 32768*3), (32768,3)).astype(np.float64)
    pos *= 1.0/65536

    _data = pos
    return _data



def test_cosmo_isol():
    """ Cosmological (non-periodic) """
    pos = get_cosmo_pos()

    for rcut, n in ndom_isol.items():
        domains = fof(pos, rcut)
        n_doms = len(np.unique(domains))
        assert(n==n_doms)
        print(rcut, n_doms, n)


def test_cosmo_periodic():
    """ Cosmological (periodic) """
    pos = get_cosmo_pos()

    for rcut, n in ndom_prd.items():
        domains = fof(pos, rcut, boxsize=1.0)
        n_doms = len(np.unique(domains))
        assert(n==n_doms)
        print(rcut, n_doms, n)
    
              
def test_shift_scale():
    """ Cube scaled and shifted """
    shift, scale = (8.82334, 4.234, -15.234), 4.8271

    pos = get_cosmo_pos()*scale + shift

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
