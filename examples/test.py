from __future__ import absolute_import, print_function
from hfof.fof import *

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
    """ Random n^3 point placement """
    from lizard.log import VerboseTimingLog
    from kdcount.cluster import fof as kdfof, dataset
    log = VerboseTimingLog()


#    for i, ngrid in [(1024, 8867), (128,1107), (256, 2217), (512, 4433)]:
    for i, ngrid in [(128,1107), (256, 2217), (512, 4433)]:
        pos = get_sim_pos(i)
        print(pos.min(), pos.max(), pos.dtype, file=log)
        rcut = (3**0.5)/ngrid
        print('rcut', rcut,file=log)

        print('Loading lib', file=log)
        cells = get_cells(pos[:5], ngrid, log)
        print('Finding cells', file=log)
        cells = get_cells(pos, ngrid, log)
        sort_idx, cellbin_data = bin_id_data(cells, log)
        sort_idx = sort_idx.astype(uint32)
#        print('Making N,3 sorted xyzw array',file=log)
#        pos = pos[sort_idx]
    
        print(cellbin_data[0].dtype,cellbin_data[1].dtype,file=log)
        

        print('3d fof', file=log)
        domains = fof3d(cellbin_data, ngrid, rcut, sort_idx, pos, log=log)
    
        print('Number of unique domains', len(unique(domains)),file=log)

        continue
        print('Trying with kdcount',file=log)
        d = dataset(pos)
        f = kdfof(d, rcut)
        print('Size', len(f.length), file=log)

    
if __name__=='__main__':
    time_kdcount()

