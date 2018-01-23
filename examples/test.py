from __future__ import absolute_import, print_function
from hfof.fof import *
import math

#from hfof.lib import qsort_pos
def get_sim_pos(i):
    
    from numpy import load, concatenate, empty, float64
    if i<=512:
        pos = load('/mainvol/peter.creasey/lizbench/data/pos_%d.npz'%i)['pos']
        return pos
    # 1024 
    spos = load('/mainvol/peter.creasey/lizbench/data/pos_512.npz')['pos']
    spos *= 0.5
 
    npos = len(spos)
    pos = empty((npos*8, 3), dtype=float64)
    pos[:npos] = spos
    pos[1*npos:2*npos] = spos+(0,0,0.5)
    pos[2*npos:3*npos] = spos+(0,0.5,0)
    pos[3*npos:4*npos] = spos+(0,0.5,0.5)
    pos[4*npos:5*npos] = spos+(0.5,0,0)
    pos[5*npos:6*npos] = spos+(0.5,0,0.5)
    pos[6*npos:7*npos] = spos+(0.5,0.5,0)
    pos[7*npos:8*npos] = spos+(0.5,0.5,0.5)
    return pos
        
def time_kdcount():
    """ Clustered n^3 points from real simulations """
    from lizard.log import VerboseTimingLog
    from kdcount.cluster import fof as kdfof, dataset
    log = VerboseTimingLog()


#    for i, ngrid in [(1024, 8867), (128,1107), (256, 2217), (512, 4433)]:
    for i, ngrid in [(128,1107), (256, 2217)]:#, (512, 4433)]:
        pos = get_sim_pos(i)
        pos_min = pos.min(axis=0)
        pos_max = pos.max(axis=0)

        box_dims = pos_max - pos_min
        max_dim = max(box_dims)
        rcut = (3**0.5)/ngrid

        cell_width = float(rcut / (3**0.5))

        inv_cell_width = float(1.0/cell_width)
        
        nval = int(math.ceil(max_dim*inv_cell_width))+2
        nval = nval | 3 # Make sure last two bits set to 1

        print('ngrid %d nval %d'%(ngrid, nval), file=log)

        print(pos_min, pos_max, pos.dtype, file=log)
        
        print('rcut', rcut,file=log)

        print('Loading lib', file=log)
        cells = get_cells(pos[:5], inv_cell_width, ngrid, log)
        print('Finding cells', file=log)
        cells = get_cells(pos, inv_cell_width, nval, log)

        print('Sorting cells', file=log)        
        sort_idx = argsort(cells)
        print('3d fof', file=log)
        domains = fof3d(cells, nval, rcut, sort_idx, pos, log=log)
    
        print('Number of unique domains', len(unique(domains)),file=log)

        print('Trying with kdcount',file=log)

        d = dataset(pos)#, boxsize=1.0)
        f = kdfof(d, rcut)
        print('Size', len(f.length), file=log)
        print('Trying with periodic kdcount',file=log)
        d = dataset(pos, boxsize=1.0)
        f = kdfof(d, rcut)
        print('Size', len(f.length), file=log)
    
if __name__=='__main__':
    time_kdcount()

