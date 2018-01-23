from __future__ import print_function, division, unicode_literals
from numpy.ctypeslib import ndpointer
from numpy.linalg import eigvalsh
import ctypes
from numpy import float64, empty, array, int32, zeros, float32, require, int64, uint32, complex128
from numpy import roll, diff, flatnonzero, uint64, cumsum, square
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
    # void find_lattice(const double *pos, const int num_pos, const int nx, int *out)
    func = _libhfof.find_lattice
    func.restype = None
    func.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_int, ndpointer(int64)]

    # Friends of Friends linking
    # int fof_link_cells(const int num_cells, const int N,const double rcut, const int64_t *restrict cell_ids, 
    #                    const int32_t *restrict cell_start_end, const uint32_t *restrict sort_idx, int32_t *restrict domains, const double *restrict xyzw)
    func = _libhfof.fof_link_cells
    func.restype = ctypes.c_int
    func.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_double, ndpointer(int64),ndpointer(int32),ndpointer(uint32), ndpointer(int32), ndpointer(float64)]


    return _libhfof

def fof3d(cell_data, ngrid, rcut, sort_idx, xyzw, log=sys.stdout):
    n = cell_data[0].shape[0]
    cell_ids = cell_data[0]
    cell_start_end = cell_data[1]
    out = empty(n, dtype=int32)
    lib = _initlib(log)
    res = lib.fof_link_cells(n, ngrid, rcut, cell_ids, cell_start_end, sort_idx, out, xyzw)

    if res<0:
        raise Exception('Error with code %d'%res)
    print('Number of domains', res, file=log)
    return out


def get_cells(pts, ncell, log=sys.stdout):
    """
    For an (N,3) array of points in [0,1), find lattice index for 
    (ncell**3,) array
    """
    lib = _initlib(log)
    p = require(pts, dtype=float64, requirements=['C']) 
    npts = p.shape[0]
    assert(p.shape ==(npts,3))
    out = empty(npts, dtype=int64)

    res = lib.find_lattice(p, npts, ncell, out)
    return out

