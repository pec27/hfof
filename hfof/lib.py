from __future__ import print_function, division, unicode_literals
from numpy.ctypeslib import ndpointer
from numpy.linalg import eigvalsh
import ctypes
from numpy import float64, empty, array, int32, zeros, float32, require, int64, uint32, complex128
from numpy import roll, diff, flatnonzero, uint64, cumsum, square, unique
from os import path
import sys

_libhfof = None

def _initlib(log):
    """ Init the library (if not already loaded) """
    global _libhfof

    if _libhfof is not None:
        return _libhfof


    name = path.join(path.dirname(path.abspath(__file__)), '../build/libhfof.so')
    if not path.exists(name):
        raise Exception('Library '+str(name)+' does not exist. Maybe you forgot to make it?')

    print('Loading libhfof - C functions for FoF calculations', name, file=log)
    _libhfof = ctypes.cdll.LoadLibrary(name)

    # Find the cell for each point
    # void find_lattice(const double *pos, const int num_pos, const double inv_cell_width, const int nx, int *out)
    func = _libhfof.find_lattice
    func.restype = None
    func.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_int, ndpointer(int64)]

    # Friends of Friends linking
    # int fof_link_cells(const int num_pos, const int N, const double b, 
    #		   const double *restrict xyz, const int64_t *restrict cells, 
    #		   const int64_t *restrict sort_idx, int32_t *restrict domains)
    func = _libhfof.fof_link_cells
    func.restype = ctypes.c_int
    func.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_double, ndpointer(float64), ndpointer(int64),ndpointer(int64), ndpointer(int32)]

    # Friends of Friends periodic linking
    # int fof_periodic(const int num_pos, const int N, const int num_orig, 
    #		 const double boxsize, const double b, 
    #		 const double *restrict xyz, const int64_t *restrict cells, 
    #		 const int64_t *restrict sort_idx, int32_t *restrict domains)
    func = _libhfof.fof_periodic
    func.restype = ctypes.c_int
    func.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double,
                     ndpointer(float64), ndpointer(int64),ndpointer(int64), ndpointer(int64), ndpointer(int32)]

    return _libhfof


def fof3d(cells, ngrid, rcut, sort_idx, xyz, log=sys.stdout):
    npos = xyz.shape[0]
    cells = require(cells, dtype=int64, requirements=['C'])
    sort_idx = require(sort_idx, dtype=int64, requirements=['C'])
    out = empty(npos, dtype=int32)
    lib = _initlib(log)
    res = lib.fof_link_cells(npos, ngrid, rcut, xyz, cells, sort_idx, out)

    if res<0:
        raise Exception('Error with code %d'%res)
    print('Number of domains {:,}'.format(res), file=log)
    return out


def fof3d_periodic(cells, ngrid, n_orig, pad_idx, rcut, sort_idx, xyz, log=sys.stdout):
    npos = xyz.shape[0]
    assert(npos>=n_orig) # cant have fewer points than non-image points
    cells = require(cells, dtype=int64, requirements=['C'])
    sort_idx = require(sort_idx, dtype=int64, requirements=['C'])
    out = empty(n_orig, dtype=int32) # only original points get domains
    lib = _initlib(log)
    res = lib.fof_periodic(npos, ngrid, n_orig, rcut, xyz, cells, sort_idx, pad_idx, out)

    if res<0:
        raise Exception('Error with code %d'%res)
    print('Number of domains {:,}'.format(res), file=log)
    return out

def get_cells(pts, inv_cell_width, ncell, log=sys.stdout):
    """
    For an (N,3) array of points in [0,1), find lattice index for 
    (ncell**3,) array
    """
    lib = _initlib(log)
    p = require(pts, dtype=float64, requirements=['C']) 
    npts = p.shape[0]
    assert(p.shape ==(npts,3))
    out = empty(npts, dtype=int64)

    res = lib.find_lattice(p, npts, inv_cell_width, ncell, out)
    return out

