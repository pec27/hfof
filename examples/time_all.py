from __future__ import absolute_import, print_function
from test import get_sim_pos
from lizard.log import VerboseTimingLog
from hfof import fof, fof_periodic
from time import time
from nbodyshop import nbodyshop_fof
from kdcount.cluster import fof as kdfof, dataset
from utils import get_sim_pos


def time_all(drt):
    """ Clustered n^3 points from real simulations """
    
    log = VerboseTimingLog()

    times = {'morton':[], 'native':[], 'random':[]}
    # Do smallest twice to average
    for i in [256]:#,512,1024]:
        for order in times.keys():
            print('Loading positions (in order=%s)'%order, file=log)
            pos = get_sim_pos(i, drt, order=order)

            rcut = 0.2/i
            time_hfof = time()
            domains = fof(pos, rcut, log=log)
            
            time_hfof = '%.3f'%(time() - time_hfof)
            
            print('Trying with nbodyshop fof',file=log)
            
            kd_setup, search, ngroups = nbodyshop_fof(pos, rcut)
            print('Times (%.3f, %.3f)'%(kd_setup, search), 'Groups {:,}'.format(ngroups), file=log)
            
            
            print('Trying with kdcount',file=log)
            
            d = dataset(pos)#, boxsize=1.0)
            time_kd = time()
            f = kdfof(d, rcut)
            time_kd = '%.3f'%(time() - time_kd)
            print('kdcount size {:,}'.format(len(f.length)), file=log)
            times[order].append((time_hfof, time_kd, '%.3f'%kd_setup, '%.3f'%search))

    print('times =', str(times), file=log)

if __name__=='__main__':
    
    from sys import argv
    drt = argv[1]
    time_all(drt)
