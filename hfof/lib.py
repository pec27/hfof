from __future__ import print_function, division, unicode_literals
from numpy.ctypeslib import ndpointer
from numpy.linalg import eigvalsh
import ctypes
from numpy import float64, empty, array, int32, zeros, float32, require, int64, uint32, complex128
from numpy import roll, diff, flatnonzero, uint64, cumsum, square, unique
from os import path
import sys

_libhfof = None

def _initlib():
    """ Init the library (if not already loaded) """
    global _libhfof

    if _libhfof is not None:
        return _libhfof


    name = path.join(path.dirname(path.abspath(__file__)), '../build/libhfof.so')
    if not path.exists(name):
        raise Exception('Library '+str(name)+' does not exist. Maybe you forgot to make it?')

    print('Loading libhfof - C functions for FoF calculations', name)
    _libhfof = ctypes.cdll.LoadLibrary(name)

    # morton indexing
    # void get_morton_idx(const double *pos, const int num_pos, const double inv_cell_width, int64_t *restrict out)
    func = _libhfof.get_morton_idx
    func.restype = None
    func.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_double, ndpointer(int64)]
    
    # minimum and maximum per cell
    # void get_min_max(const double *pos, const int num_pos, double *restrict out)
    func = _libhfof.get_min_max
    func.restype = None
    func.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ndpointer(ctypes.c_double)]
    
    # Find the cell for each point
    # void find_lattice(const double *pos, const int num_pos, 
    #                   const double inv_cell_width, const int N, const int M, int64_t *out)
    func = _libhfof.find_lattice
    func.restype = None
    func.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                     ctypes.c_int, ctypes.c_int, ndpointer(int64)]

    # Find the block+cell for each point
    # void blocks_cells(const double min_x, const double min_y, const double min_z, 
    #		  const double *pos, const int num_pos, 
    #		  const double inv_cell_width, const int ny, const int nx, 
    #		  int64_t *out)
    func = _libhfof.blocks_cells
    func.restype = None
    func.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double,
                     ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                     ctypes.c_int, ctypes.c_int, ndpointer(int64)]

    # Friends of Friends linking periodic (on 4x4x4 cells)
    # see src/fof64.c
    func = _libhfof.fof64
    func.restype = ctypes.c_int
    func.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_int64,ctypes.c_int,ctypes.c_double, 
                     ndpointer(float64), ndpointer(int64),ndpointer(int64), ndpointer(int64), 
                     ndpointer(int32), ctypes.c_double]

    # Friends of Friends periodic linking
    # see src/fof.c
    func = _libhfof.fof_periodic
    func.restype = ctypes.c_int
    func.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_int64,ctypes.c_int,ctypes.c_double,
                     ndpointer(float64), ndpointer(int64),ndpointer(int64), ndpointer(int64), ndpointer(int32)]

    # Periodic image insertion
    # int pad_box(const double inv_boxsize, const double r_pad, const int num_pos, 
    #          const double *restrict pos, double *restrict periodic_pos)
    # 	      int64_t *restrict pad_idx, const int max_images)

    func = _libhfof.pad_box
    func.restype = ctypes.c_int
    func.argtypes = [ctypes.c_double, ctypes.c_double,ctypes.c_int,
                     ndpointer(float64), ndpointer(float64), ndpointer(int64), ctypes.c_int]
    
    return _libhfof

def fof3d_periodic(cells, N, M, n_orig, pad_idx, rcut, sort_idx, xyz, log=None):
    npos = xyz.shape[0]
    assert(npos>=n_orig) # cant have fewer points than non-image points
    cells = require(cells, dtype=int64, requirements=['C'])
    sort_idx = require(sort_idx, dtype=int64, requirements=['C'])
    out = empty(n_orig, dtype=int32) # only original points get domains
    lib = _initlib()
    res = lib.fof_periodic(npos, int(N), int(M), int(n_orig), rcut, xyz, cells, sort_idx, pad_idx, out)

    if res<0:
        raise Exception('Error with code %d'%res)
    if log is not None:
        print('Number of domains {:,}'.format(res), file=log)
    return out

def fof_periodic64(cells, N, M, rcut, sort_idx, xyz, periodic_pad_idx=None, log=None):
    npos = xyz.shape[0]
    if periodic_pad_idx is None:
        n_orig = npos
        pad_idx = empty(0, dtype=int64) # Dummy, no images
    else:
        n_orig = npos - len(periodic_pad_idx)
        pad_idx = require(periodic_pad_idx, dtype=int64, requirements=['C'])

    cells = require(cells, dtype=int64, requirements=['C'])
    sort_idx = require(sort_idx, dtype=int64, requirements=['C'])
    out = empty(n_orig, dtype=int32) # only original points get domains
    desired_load=0.6
    lib = _initlib()
    res = lib.fof64(npos, int(N), int(M), int(n_orig), rcut, xyz, cells, sort_idx, pad_idx, out, desired_load)

    if res<0:
        raise Exception('Error with code %d'%res)
    if log is not None:
        print('Number of domains {:,}'.format(res), file=log)
    return out

def get_cells(pts, inv_cell_width, Ny, Nx, log=sys.stdout):
    """
    For an (N,3) array of points in [0,1), find lattice index for 
    (ncell**3,) array
    """
    lib = _initlib()
    p = require(pts, dtype=float64, requirements=['C']) 
    npts = p.shape[0]
    assert(p.shape ==(npts,3))
    out = empty(npts, dtype=int64)

    res = lib.find_lattice(p, npts, inv_cell_width, Ny, Nx, out)
    return out


def get_blocks_cells(pts, inv_cell_width, Ny, Nx, pts_min, log=sys.stdout):
    """
    For an (N,3) array of points in [0,1), find index 
    idx = block<<6 | cell
    where 
    ix,iy,iz = (int)(pts[x,y,z]*inv_cell_width)
    block := (ix>>2)*Nx + (iy>>2)*Ny + (iz>>2) 
    cell := (ix&3)<<4 | (iy&3)<<2 | (iz&0x3) # (i.e. 0-63)
    """
    lib = _initlib()
    p = require(pts, dtype=float64, requirements=['C']) 
    npts = p.shape[0]
    assert(p.shape ==(npts,3))

    out = empty(npts, dtype=int64)
    min_x, min_y, min_z = pts_min.astype(float64)

    res = lib.blocks_cells(min_x, min_y, min_z, p, npts, inv_cell_width, 
                           Ny, Nx, out)
    return out

def minmax(pts):
    lib = _initlib()
    p = require(pts, dtype=float64, requirements=['C']) 
    out = empty(6, dtype=float64)
    res = lib.get_min_max(p, len(p), out)
    return out[:3], out[3:]

def morton_idx(pts):
    """ Find the morton index of the points (on an 8192^3 grid) """
    lib = _initlib()
    p = require(pts, dtype=float64, requirements=['C']) 
    inv_cell_width = 1.0/8192
    npts = len(p)
    out = empty(npts, dtype=int64)
    lib.get_morton_idx(p, npts, inv_cell_width, out)
    return out
    
def pad_scaled_cube(pts, boxsize, r_pad, log=None):
    """

    For a set of points, find images in [0,boxsize)^3, scale this [0,1)^3 and
    then add the right-hand repeats in [1,r_pad/boxsize) for each dimension
    (incl. repeats of repeats). Return arrays of new positions (scaled and 
    padded) along with the corresponding indices in the original array. The 
    scaled points are always the first N pts in the new array

    pts     - (N,3) array
    boxsize - size for periodicity
    r_pad   - edge to pad onto [boxsize, boxsize+r_pad)
    
    returns pad_idx, pos
       pad_idx - (N_new) indices of the positions used to pad
       new_pos - (N+N_new, ndim) array of scaled pos + padded pos [0, 1+r_pad/boxsize]

    """
    lib = _initlib()
    p = require(pts, dtype=float64, requirements=['C']) 
    inv_boxsize = 1.0/boxsize
    N = len(p)

    fudge = 2.0 
    guess_images = int(fudge*N*3.0*r_pad/boxsize) + 1
    
    new_pos = empty((N+guess_images,3), dtype=float64)
    pad_idx = empty((guess_images,), dtype=int64)
    
    new_ims = lib.pad_box(inv_boxsize, r_pad, N, p, new_pos, pad_idx, guess_images)
    while new_ims<0:
        guess_images *= 2
        if log is not None:
            print('Too many images! Doubling size...', file=log)

        new_pos = empty((N+guess_images,3), dtype=float64)
        pad_idx = empty((guess_images,), dtype=int64)
        new_ims = lib.pad_box(inv_boxsize, r_pad, N, p, new_pos, pad_idx, guess_images)
        
    if log is not None:
        print('{:,} images (c.f. guessed={:,})'.format(new_ims, guess_images), file=log)

    # Arrays contain only used vals
    new_pos = new_pos[:N+new_ims] # originals+images
    pad_idx = pad_idx[:new_ims] # images
    return pad_idx, new_pos



