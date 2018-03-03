from os import path
from numpy import load, concatenate, empty, float64
import numpy as np
from numpy.random import shuffle
from hfof.lib import morton_idx

def get_sim_pos(i, drt = '/mainvol/peter.creasey/lizbench/data/', order='native'):
    
    pos = 0
    if i<=512:
        pos = load(path.join(drt,'pos_%d.npz'%i))['pos']
    else:
        # 1024 
        spos = load(path.join(drt,'pos_512.npz'))['pos']
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

    if order=='morton':
        idx = morton_idx(pos)
        idx = np.argsort(idx)
        pos = pos[idx]
    elif order=='random':
        shuffle(pos)
    
    return pos
