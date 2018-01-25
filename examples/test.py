from __future__ import absolute_import, print_function
from hfof.fof import *



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

    print('Loading lib', file=log)
    cells = get_cells([(0.5, 0.2, 0.1)], 0.333, 3, log)

#    for i, ngrid in [(1024, 8867), (128,1107), (256, 2217), (512, 4433)]:
    for i, ngrid in [(128,1107), (256, 2217)]:#, (512, 4433)]:
        pos = get_sim_pos(i)

        rcut = (3**0.5)/ngrid

        domains = fof(pos, rcut, log=log)

        print('Number of unique domains {:,}'.format(len(unique(domains))),file=log)
        continue

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

